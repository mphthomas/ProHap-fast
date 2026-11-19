#!/usr/bin/env python3
"""Simulate a phased, WGS-like multi-sample VCF over a block of MANE Select transcripts, for benchmarks and tests.

    python3 bench/simulate.py --refs ~/prohap-refs/ensembl.v113 --chrom 22 --first ENST00000402142 \\
        --transcripts 60 --samples 26000 --out DIR

writes into DIR:

    chr22.vcf.gz       the VCF (gzip level 1); with --chunks K, chr22.1.vcf.gz ... chr22.K.vcf.gz instead, so the
                       workflow's glob is chr{chr}.*vcf.gz
    samples.tsv        sample metadata in the layout of igsr_samples.tsv: name, sex, 26 1000 Genomes populations
    transcripts.csv    the block's transcripts (chromosome,transcriptID); rows of other chromosomes already in the
                       file are kept, so several runs into one directory build one dataset
    chr22.json         the arguments and what was simulated

The output depends only on the arguments and the references: the same command writes the same bytes. Each part has
its own random stream (sample metadata from --seed alone, sites and genotypes from --seed and the chromosome, IDs
from a third stream), so runs for different chromosomes share one samples.tsv, and --ids changes only the ID column.

The block is --transcripts consecutive MANE Select transcripts of the chromosome (in order of start, only those with
an annotated start codon, as prohap.py keeps them; 0 for all to the end), starting at --first or after skipping
--skip. REF bases come from
the Ensembl cDNA mapped back through the exons (complemented on the - strand), so REF matches the reference that
prohap.py checks. The model (docs/IMPLEMENTATION_NOTES.md, section 5):

- sites per exonic bp: 0.3 * (N / 500,000) ** 0.43 for N samples, without replacement over the block's exonic bases;
- common sites, one per 300 exonic bp: allele frequency log-uniform in 0.01-0.5, carried by founder haplotypes. Each
  50 kb window has 40 founders with Dirichlet(0.5) weights; every haplotype copy descends from one founder per window,
  and a site is carried by a random set of founders whose weights add up to its frequency. This gives shared
  haplotypes, as in real cohorts;
- the other sites are rare: allele count from a discrete power law (exponent 1.9) up to 1% of 2N, carried by random
  copies;
- 6% deletions and 6% insertions (1-10 bp, geometric), 6% with a second ALT allele (rare, on other copies); a
  deletion may run past the end of its exon (the intronic REF bases are random: nothing reads them);
- 10% of deletions have an SNV inside the deleted bases on some of the same copies, so prohap.py finds conflicting
  variants and logs removed samples;
- IDs: `dot` (every ID is `.`; prohap.py names them by line number), `rs` (unique rs numbers), or `ambiguous` (30%
  tokens pandas may retype: 2e5 123 7 0012 NA nan True False 1.50, the rest `.`);
- chromosome X outside the pseudoautosomal regions: males carry one copy, written haploid (`1`) or mirrored (`1|1`);
- FORMAT is GT only; INFO is AC, AN and AF (6 significant digits), as in 1000 Genomes.
"""
import argparse
import gzip
import hashlib
import io
import json
import os
import sys
import time

import numpy as np

# 26 populations of the 1000 Genomes Project (names as in igsr_samples.tsv)
POPULATIONS = [
    ('ACB', 'African Caribbean', 'AFR'), ('ASW', 'African Ancestry SW', 'AFR'), ('BEB', 'Bengali', 'SAS'),
    ('CDX', 'Dai Chinese', 'EAS'), ('CEU', 'CEPH', 'EUR'), ('CHB', 'Han Chinese', 'EAS'),
    ('CHS', 'Southern Han Chinese', 'EAS'), ('CLM', 'Colombian', 'AMR'), ('ESN', 'Esan', 'AFR'),
    ('FIN', 'Finnish', 'EUR'), ('GBR', 'British', 'EUR'), ('GIH', 'Gujarati', 'SAS'),
    ('GWD', 'Gambian Mandinka', 'AFR'), ('IBS', 'Iberian', 'EUR'), ('ITU', 'Telugu', 'SAS'),
    ('JPT', 'Japanese', 'EAS'), ('KHV', 'Kinh Vietnamese', 'EAS'), ('LWK', 'Luhya', 'AFR'), ('MSL', 'Mende', 'AFR'),
    ('MXL', 'Mexican Ancestry', 'AMR'), ('PEL', 'Peruvian', 'AMR'), ('PJL', 'Punjabi', 'SAS'),
    ('PUR', 'Puerto Rican', 'AMR'), ('STU', 'Tamil', 'SAS'), ('TSI', 'Toscani', 'EUR'), ('YRI', 'Yoruba', 'AFR'),
]
SUPERPOPULATIONS = {'AFR': 'African Ancestry', 'AMR': 'American Ancestry', 'EAS': 'East Asian Ancestry',
                    'EUR': 'European Ancestry', 'SAS': 'South Asian Ancestry'}

# the model
DENSITY_AT_500K = 0.3       # sites per exonic bp at 500,000 samples ...
DENSITY_EXPONENT = 0.43     # ... scaled by (N / 500,000) ** 0.43
COMMON_EVERY = 300          # one common site per this many exonic bp
COMMON_AF = (0.01, 0.5)     # log-uniform
FOUNDERS = 40
FOUNDER_WINDOW = 50_000
DIRICHLET_ALPHA = 0.5
RARE_EXPONENT = 1.9
RARE_MAX_AF = 0.01
P_DELETION = 0.06
P_INSERTION = 0.06
P_SECOND_ALT = 0.06
INDEL_P = 0.5               # geometric length distribution, 1..INDEL_MAX
INDEL_MAX = 10
P_NESTED_SNV = 0.1          # deletions with an SNV inside on some of the same copies
AMBIGUOUS_FRACTION = 0.3
X_PAR1_TO = 2781479         # the workflow's defaults (x_par1_to, x_par2_from)
X_PAR2_FROM = 155701383

BASES = 'ACGT'
COMPLEMENT = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
CHROM_CODE = {**{str(i): i for i in range(1, 23)}, 'X': 23, 'Y': 24, 'MT': 25}

# stream identifiers for np.random.SeedSequence([seed, stream, chromosome])
STREAM_SAMPLES, STREAM_SITES, STREAM_IDS = 0, 1, 2


def log(*args):
    print('[simulate]', *args, file=sys.stderr, flush=True)


def default_refs():
    root = os.environ.get('PROHAP_REFS', '')
    found = sorted(d for d in (os.listdir(root) if root and os.path.isdir(root) else []) if d.startswith('ensembl.v'))
    return os.path.join(root, found[0]) if len(found) == 1 else None


def read_release(refs):
    with open(os.path.join(refs, 'release.txt')) as fh:
        return fh.read().strip()


# ---------------------------------------------------------------- the block

def block_transcripts(refs, release, chrom, first, skip, count):
    """The block: [(transcript feature, exons sorted by start)], consecutive MANE Select transcripts with a start
    codon, in order of start (as prohap.py sorts them)."""
    import gffutils
    db = gffutils.FeatureDB(os.path.join(refs, 'Homo_sapiens.GRCh38.%s.chr_patch_hapl_scaff_chr%s.db'
                                         % (release, chrom)))
    with open(os.path.join(refs, 'transcripts_reference_%s.csv' % release)) as fh:
        ids = [line.rstrip('\n').split(',')[1] for line in fh if line.split(',', 1)[0] == chrom]
    features = []
    for tid in ids:
        feature = db[tid]
        if next(iter(db.children(feature, featuretype='start_codon')), None) is not None:
            features.append(feature)
    features.sort(key=lambda f: f.start)
    if first:
        starts = [i for i, f in enumerate(features) if f.id == first]
        if not starts:
            sys.exit('simulate: %s is not a MANE Select transcript with a start codon on chromosome %s'
                     % (first, chrom))
        skip = starts[0]
    chosen = features[skip:skip + count] if count else features[skip:]
    if len(chosen) < max(count, 1):
        sys.exit('simulate: chromosome %s has only %d MANE Select transcripts with a start codon after skipping %d'
                 % (chrom, len(chosen), skip))
    return [(f, list(db.children(f, featuretype='exon', order_by='start'))) for f in chosen]


def read_cdnas(path, wanted):
    """{transcript ID without version: sequence} for the wanted IDs."""
    found = {}
    current = None
    parts = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if current is not None:
                    found[current] = ''.join(parts)
                tid = line[1:].split(None, 1)[0].split('.')[0]
                current = tid if tid in wanted else None
                parts = []
            elif current is not None:
                parts.append(line.strip())
    if current is not None:
        found[current] = ''.join(parts)
    return found


def exonic_bases(block, cdnas):
    """{genomic position: base on the + strand} for every exonic position of the block whose transcript's cDNA is
    the exons spliced together; and the transcripts whose cDNA is not."""
    base = {}
    unusable = []
    for feature, exons in block:
        cdna = cdnas.get(feature.id.split('.')[0])
        if cdna is None or len(cdna) != sum(e.end - e.start + 1 for e in exons):
            unusable.append(feature.id)
            continue
        plus = cdna if feature.strand == '+' else cdna.translate(COMPLEMENT)[::-1]
        offset = 0
        for e in exons:
            for i in range(e.end - e.start + 1):
                base.setdefault(e.start + i, plus[offset + i].upper())
            offset += e.end - e.start + 1
    return base, unusable


# ---------------------------------------------------------------- samples

def sample_names(n, prefix='S'):
    width = max(6, len(str(n)))
    return ['%s%0*d' % (prefix, width, i + 1) for i in range(n)]


def make_samples(n, seed, unlisted):
    """VCF sample names, sexes (True = male), and the metadata rows (unlisted VCF samples left out, as many
    metadata-only samples added)."""
    rng = np.random.default_rng(np.random.SeedSequence([seed, STREAM_SAMPLES, n, unlisted]))
    names = sample_names(n)
    male = rng.random(n) < 0.5
    pop = rng.integers(0, len(POPULATIONS), n)
    left_out = set(rng.choice(n, size=unlisted, replace=False).tolist()) if unlisted else set()
    rows = [(names[i], male[i], POPULATIONS[pop[i]]) for i in range(n) if i not in left_out]
    extra_male = rng.random(unlisted) < 0.5
    extra_pop = rng.integers(0, len(POPULATIONS), unlisted)
    rows += [(name, extra_male[i], POPULATIONS[extra_pop[i]]) for i, name in enumerate(sample_names(unlisted, 'M'))]
    return names, male, rows


def write_samples(path, rows):
    text = 'Sample name\tSex\tPopulation code\tPopulation name\tSuperpopulation code\tSuperpopulation name\n' + ''.join(
        '%s\t%s\t%s\t%s\t%s\t%s\n' % (name, 'male' if m else 'female', p[0], p[1], p[2], SUPERPOPULATIONS[p[2]])
        for name, m, p in rows)
    if os.path.exists(path):
        with open(path) as fh:
            if fh.read() != text:
                sys.exit('simulate: %s exists with other samples; use the same --samples, --seed and --unlisted for '
                         'every chromosome of a dataset, or another --out' % path)
        return
    with open(path, 'w') as fh:
        fh.write(text)


# ---------------------------------------------------------------- sites and genotypes

def indel_length(rng):
    return min(int(rng.geometric(INDEL_P)), INDEL_MAX)


def random_bases(rng, n):
    return ''.join(BASES[i] for i in rng.integers(0, 4, n))


def other_base(rng, ref, avoid=''):
    return rng.choice([b for b in BASES if b != ref and b not in avoid])


def rare_count_sampler(n_copies):
    top = max(1, int(RARE_MAX_AF * n_copies))
    counts = np.arange(1, top + 1)
    p = counts.astype(float) ** -RARE_EXPONENT
    return counts, p / p.sum()


def simulate_sites(rng, base, n, male, chrom):
    """[(pos, ref, [alt, ...], [carrier copies per alt])] sorted by position. Copy 2i is sample i's first haplotype
    (left of the |), 2i + 1 its second."""
    positions = np.array(sorted(base), dtype=np.int64)
    n_bp = len(positions)
    n_sites = min(n_bp, int(round(DENSITY_AT_500K * (n / 500000) ** DENSITY_EXPONENT * n_bp)))
    n_common = min(n_sites, int(round(n_bp / COMMON_EVERY)))
    site_pos = np.sort(rng.choice(positions, size=n_sites, replace=False))
    is_common = np.zeros(n_sites, dtype=bool)
    is_common[rng.choice(n_sites, size=n_common, replace=False)] = True
    kind = rng.random(n_sites)
    taken = set(site_pos.tolist())

    n_copies = 2 * n
    eligible_all = np.ones(n_copies, dtype=bool)
    eligible_x = eligible_all.copy()
    eligible_x[2 * np.flatnonzero(male) + 1] = False        # males' second copy outside the PARs

    def eligible(pos):
        if chrom == 'X' and X_PAR1_TO < pos < X_PAR2_FROM:
            return eligible_x
        return eligible_all

    rare_counts, rare_p = rare_count_sampler(n_copies)

    # founders: per 50 kb window, the founder of every copy
    start = int(site_pos[0]) if n_sites else 0
    founders = {}

    def founder_of(pos):
        w = (pos - start) // FOUNDER_WINDOW
        if w not in founders:
            weights = rng.dirichlet(np.full(FOUNDERS, DIRICHLET_ALPHA))
            founders[w] = (rng.choice(FOUNDERS, size=n_copies, p=weights).astype(np.uint8), weights)
        return founders[w]

    def common_carriers(pos, ok):
        assigned, weights = founder_of(pos)
        af = float(np.exp(rng.uniform(np.log(COMMON_AF[0]), np.log(COMMON_AF[1]))))
        chosen = []
        total = 0.0
        for f in rng.permutation(FOUNDERS):
            if total >= af:
                break
            chosen.append(f)
            total += weights[f]
        carriers = np.flatnonzero(np.isin(assigned, chosen) & ok)
        return carriers if len(carriers) else rare_carriers(ok)     # founders too rare in a small cohort

    def rare_carriers(ok, exclude=None):
        pool = np.flatnonzero(ok if exclude is None else ok & ~exclude)
        c = min(int(rng.choice(rare_counts, p=rare_p)), len(pool))
        return np.sort(rng.choice(pool, size=c, replace=False))

    sites = []
    for i in range(n_sites):
        pos = int(site_pos[i])
        ref = base[pos]
        ok = eligible(pos)
        carriers = common_carriers(pos, ok) if is_common[i] else rare_carriers(ok)
        if kind[i] < P_DELETION:
            length = indel_length(rng)
            ref = ''.join(base.get(pos + k, None) or random_bases(rng, 1) for k in range(length + 1))
            sites.append((pos, ref, [ref[0]], [carriers]))
            if rng.random() < P_NESTED_SNV:
                inside = [q for q in range(pos, pos + length + 1) if q in base and (q == pos or q not in taken)]
                if inside and len(carriers):
                    q = int(rng.choice(inside))
                    taken.add(q)
                    shared = rng.choice(carriers, size=max(1, len(carriers) // 2), replace=False)
                    extra = rng.choice(np.flatnonzero(eligible(q)), size=int(rng.integers(0, 3)), replace=False)
                    nested = np.unique(np.concatenate([shared, extra]))
                    sites.append((q, base[q], [other_base(rng, base[q])], [nested]))
        elif kind[i] < P_DELETION + P_INSERTION:
            sites.append((pos, ref, [ref + random_bases(rng, indel_length(rng))], [carriers]))
        elif kind[i] < P_DELETION + P_INSERTION + P_SECOND_ALT:
            alt1 = other_base(rng, ref)
            alt2 = other_base(rng, ref, alt1)
            taken_mask = np.zeros(n_copies, dtype=bool)
            taken_mask[carriers] = True
            sites.append((pos, ref, [alt1, alt2], [carriers, rare_carriers(ok, taken_mask)]))
        else:
            sites.append((pos, ref, [other_base(rng, ref)], [carriers]))
    sites.sort(key=lambda s: s[0])      # stable: a nested SNV at a deletion's position comes after it
    return sites


# ---------------------------------------------------------------- IDs

def ambiguous_token(rng):
    kind = int(rng.integers(0, 9))
    if kind == 0:
        return '%de%d' % (rng.integers(1, 10), rng.integers(1, 10))    # 2e5: a float
    if kind == 1:
        return str(rng.integers(10, 1000000))                           # 123: an integer
    if kind == 2:
        return str(rng.integers(0, 10))                                 # 7
    if kind == 3:
        return '00' + str(rng.integers(1, 100))                         # 0012: an integer, 12
    if kind == 4:
        return 'NA'                                                     # NaN
    if kind == 5:
        return 'nan'
    if kind == 6:
        return 'True'                                                   # a boolean
    if kind == 7:
        return 'False'
    return '%d.%d0' % (rng.integers(0, 10), rng.integers(1, 10))        # 1.50: a float, 1.5


def make_ids(rng, mode, n_records):
    if mode == 'dot':
        return ['.'] * n_records
    if mode == 'rs':
        numbers = set()
        while len(numbers) < n_records:
            numbers.update(rng.integers(1, 1_500_000_000, n_records - len(numbers)).tolist())
        return ['rs%d' % x for x in rng.permutation(sorted(numbers))]
    ambiguous = rng.random(n_records) < AMBIGUOUS_FRACTION
    return [ambiguous_token(rng) if a else '.' for a in ambiguous]


# ---------------------------------------------------------------- writing

def format_af(x):
    return '%.6g' % x


def vcf_header(chrom, names, command):
    lines = [
        '##fileformat=VCFv4.2',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##source=ProHap-fast bench/simulate.py',
        '##simulateCommand=' + command,
        '##contig=<ID=%s>' % chrom,
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">',
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency, AC / AN">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype">',
        '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + names),
    ]
    return ('\n'.join(lines) + '\n').encode()


def write_vcfs(out, chrom, sites, ids, names, male, chunks, haploid_males, command):
    n = len(names)
    template = np.frombuffer(b'0|0\t' * n, dtype=np.uint8).copy()
    template[-1] = ord('\n')
    male_rows = np.flatnonzero(male)
    keep_haploid = np.ones(4 * n, dtype=bool)          # drop '|0' of males: '1|0\t' -> '1\t'
    keep_haploid[4 * male_rows + 1] = False
    keep_haploid[4 * male_rows + 2] = False
    an_x = 2 * n - len(male_rows)
    header = vcf_header(chrom, names, command)

    for name in os.listdir(out):
        if name == 'chr%s.vcf.gz' % chrom or (name.startswith('chr%s.' % chrom) and name.endswith('.vcf.gz')):
            os.remove(os.path.join(out, name))
    bounds = [len(sites) * k // chunks for k in range(chunks + 1)]
    paths = []
    for k in range(chunks):
        path = os.path.join(out, 'chr%s.vcf.gz' % chrom if chunks == 1 else 'chr%s.%d.vcf.gz' % (chrom, k + 1))
        paths.append(path)
        with open(path, 'wb') as raw, gzip.GzipFile(filename='', mode='wb', compresslevel=1, fileobj=raw,
                                                    mtime=0) as gz, io.BufferedWriter(gz, 1 << 22) as fh:
            fh.write(header)
            for i in range(bounds[k], bounds[k + 1]):
                pos, ref, alts, carriers = sites[i]
                non_par = chrom == 'X' and X_PAR1_TO < pos < X_PAR2_FROM
                gt = template.copy()
                for allele, copies in enumerate(carriers, start=1):
                    gt[4 * (copies // 2) + 2 * (copies % 2)] = ord(str(allele))
                # AC and AN as a VCF writer counts them from the genotypes written
                if non_par and haploid_males:
                    an, ac = an_x, [len(c) for c in carriers]
                elif non_par:       # mirrored: a male carrier counts twice
                    an, ac = 2 * n, [len(c) + int(male[c // 2].sum()) for c in carriers]
                else:
                    an, ac = 2 * n, [len(c) for c in carriers]
                info = 'AC=%s;AN=%d;AF=%s' % (','.join(map(str, ac)), an, ','.join(format_af(a / an) for a in ac))
                fh.write(('%s\t%d\t%s\t%s\t%s\t.\tPASS\t%s\tGT\t' % (chrom, pos, ids[i], ref, ','.join(alts),
                                                                     info)).encode())
                if non_par and haploid_males:
                    fh.write(gt[keep_haploid].tobytes())
                else:
                    if non_par:     # mirrored: the male's second copy repeats the first
                        gt[4 * male_rows + 2] = gt[4 * male_rows]
                    fh.write(gt.tobytes())
    return paths


def update_transcripts(path, chrom, block):
    rows = []
    if os.path.exists(path):
        with open(path) as fh:
            rows = [line.rstrip('\n') for line in fh][1:]
    rows = [r for r in rows if r.split(',', 1)[0] != chrom] + ['%s,%s' % (chrom, f.id) for f, _ in block]
    with open(path, 'w') as fh:
        fh.write('chromosome,transcriptID\n' + ''.join(r + '\n' for r in rows))


# ---------------------------------------------------------------- main

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--refs', default=default_refs(),
                        help='references directory ensembl.v<release> (default: the one in $PROHAP_REFS)')
    parser.add_argument('--chrom', required=True, help='chromosome, e.g. 22 or X')
    where = parser.add_mutually_exclusive_group()
    where.add_argument('--first', help='first transcript of the block (default: the first on the chromosome)')
    where.add_argument('--skip', type=int, default=0, help='transcripts to skip before the block (default 0)')
    parser.add_argument('--transcripts', type=int, default=60,
                        help='transcripts in the block (default 60; 0: to the end of the chromosome)')
    parser.add_argument('--samples', type=int, required=True, help='samples in the VCF')
    parser.add_argument('--seed', type=int, default=0, help='random seed (default 0)')
    parser.add_argument('--ids', choices=['dot', 'rs', 'ambiguous'], default='dot', help='ID column (default dot)')
    parser.add_argument('--males', choices=['haploid', 'mirrored'], default='haploid',
                        help='genotypes of males on X outside the PARs (default haploid)')
    parser.add_argument('--unlisted', type=int, default=0,
                        help='VCF samples left out of samples.tsv, and as many samples.tsv rows not in the VCF')
    parser.add_argument('--chunks', type=int, default=1, help='split the VCF into this many files (default 1)')
    parser.add_argument('--out', required=True, help='output directory')
    args = parser.parse_args(argv)

    if not args.refs:
        parser.error('--refs is required ($PROHAP_REFS does not hold exactly one ensembl.v<release>)')
    if args.chrom not in CHROM_CODE:
        parser.error('--chrom must be 1-22, X, Y or MT')
    if args.samples < 1 or args.transcripts < 0 or args.chunks < 1 or not 0 <= args.unlisted < args.samples:
        parser.error('--samples and --chunks must be positive, --transcripts not negative, and --unlisted below '
                     '--samples')
    refs = os.path.abspath(os.path.expanduser(args.refs))
    release = read_release(refs)
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()

    block = block_transcripts(refs, release, args.chrom, args.first, args.skip, args.transcripts)
    cdnas = read_cdnas(os.path.join(refs, 'total_cdnas_%s.fa' % release), {f.id.split('.')[0] for f, _ in block})
    base, unusable = exonic_bases(block, cdnas)
    if unusable:
        log('no matching cDNA for %s: their exons get no sites unless another transcript covers them'
            % ','.join(unusable))
    names, male, metadata = make_samples(args.samples, args.seed, args.unlisted)
    write_samples(os.path.join(args.out, 'samples.tsv'), metadata)

    code = CHROM_CODE[args.chrom]
    sites = simulate_sites(np.random.default_rng(np.random.SeedSequence([args.seed, STREAM_SITES, code])),
                           base, args.samples, male, args.chrom)
    ids = make_ids(np.random.default_rng(np.random.SeedSequence([args.seed, STREAM_IDS, code])), args.ids, len(sites))
    command = ' '.join(['simulate.py', '--chrom', args.chrom, '--first', block[0][0].id, '--transcripts',
                        str(args.transcripts), '--samples', str(args.samples), '--seed', str(args.seed), '--ids',
                        args.ids, '--males', args.males, '--unlisted', str(args.unlisted), '--chunks', str(args.chunks),
                        '(Ensembl %s)' % release])
    paths = write_vcfs(args.out, args.chrom, sites, ids, names, male, args.chunks, args.males == 'haploid', command)
    update_transcripts(os.path.join(args.out, 'transcripts.csv'), args.chrom, block)

    with open(os.path.abspath(__file__), 'rb') as fh:
        source = hashlib.sha256(fh.read()).hexdigest()
    summary = {
        'command': command, 'simulate_sha256': source, 'ensembl_release': release,
        'block': [f.id for f, _ in block], 'block_start': block[0][0].start, 'block_end': max(f.end for f, _ in block),
        'exonic_bp': len(base), 'records': len(sites),
        'af_at_least_0.01': sum(1 for s in sites if len(s[3][0]) >= COMMON_AF[0] * 2 * args.samples),
        'multiallelic': sum(1 for s in sites if len(s[2]) > 1),
        'deletions': sum(1 for s in sites if len(s[1]) > 1), 'insertions': sum(1 for s in sites if len(s[2][0]) > 1),
        'vcf': [os.path.basename(p) for p in paths], 'seconds': round(time.time() - t0, 1),
    }
    with open(os.path.join(args.out, 'chr%s.json' % args.chrom), 'w') as fh:
        json.dump(summary, fh, indent=1)
        fh.write('\n')
    log('chr%s: %d transcripts, %d exonic bp, %d records, %d samples -> %s (%.1f s)'
        % (args.chrom, len(block), len(base), len(sites), args.samples, ', '.join(summary['vcf']),
           summary['seconds']))


if __name__ == '__main__':
    main()
