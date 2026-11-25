#!/usr/bin/env python3
"""A naive, independent check of prohap.py -fix_bugs 1: the protein of every haplotype copy, computed copy by copy with
plain string edits and Biopython, against the protein the fast code gives that copy.

    python3 tests/fix_oracle.py [--work DIR] [--keep] [--chromosomes 21,22] DATASET...   # 1kg, mix, ttn, block, all

For each dataset and chromosome: the reference's filter (legacy/src/vcf_filter_fix.py) writes the filtered VCF, and
src/prohap.py -fix_bugs 1 runs on it with every haplotype kept (thresholds 0, ignore_UTR 1, skip_start_lost 1,
require_start 1). The oracle then reads the filtered VCF itself and, for each transcript and each haplotype copy
(both copies of every sample in the metadata; one for males outside the pseudo-autosomal regions of X), works out:

1. the records the copy carries that overlap an exon of the transcript;
2. conflicts (IMPLEMENTATION_NOTES section 9.2): each record reduced to the bases it changes (shared prefix, then
   shared suffix of REF and ALT removed); a copy whose changes overlap, with an insertion strictly inside a change
   or two insertions at one point, is expected in the removed-samples log;
3. the liftover check (9.3): a record whose trimmed ALT is in the reference cDNA where its trimmed REF is not is
   left out;
4. the haplotype's cDNA: the transcript's exonic bases in genome order, each trimmed change applied as a string
   edit (bases removed, ALT put in their place, inserted bases put between their two bases), reverse complemented
   on the - strand;
5. the start codon (9.5): lost when a change touches it (an insertion strictly inside it), shifted by the length
   changes before it; when lost, the first ATG of the haplotype's cDNA from where the lost start was, in any frame;
   none: no protein;
6. the protein: Biopython's translation from the start to the first stop.

A copy whose protein equals the reference protein (the reference cDNA translated the same way) is expected to have
no row; any other copy a row whose FASTA sequence is its protein. Copies with a record that the oracle does not model
are left out and counted: records across an exon boundary or an insertion at one (the reference's truncation at
splice sites decides those), and REF alleles not in the cDNA (the haplotype is skipped, case 12).

Differences are reported by kind. Known kinds that the fixes do not address are named in the report:
`shared_id` (the copy's row merges haplotypes with equal VCF ID strings, IMPLEMENTATION_NOTES 7.3).
"""
import argparse
import collections
import gzip
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REFS = Path(os.environ.get('PROHAP_REFS') or '/opt/prohap/refs').expanduser() / 'ensembl.v113'
PAR1_TO, PAR2_FROM = 2781479, 155701383
COMPLEMENT = str.maketrans('ACGTN', 'TGCAN')
NA_TOKENS = {'NA', 'nan', 'NaN', 'N/A', 'n/a', 'NULL', 'null', 'None', '<NA>', '#N/A', '-nan', '-NaN', ''}


def rc(seq):
    return seq.translate(COMPLEMENT)[::-1]


def translate(seq):
    """Biopython, whole codons only, up to the first stop (excluded)."""
    from Bio.Seq import Seq
    seq = seq[:len(seq) // 3 * 3]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        protein = str(Seq(seq).translate())
    return protein.split('*', 1)[0]


# ---------------------------------------------------------------------------------------------- datasets

def datasets(work):
    """name -> list of (chromosome, VCF glob, AF threshold, samples.tsv, transcripts.csv)."""
    sim = ROOT / 'tests' / 'fixtures' / 'sim' / 'cache'
    block = ROOT / 'bench' / 'data' / 'chr22-ENST00000402142-60-n1000-dot-seed0'

    def onekg():
        data = Path(work) / 'input-1kg'
        if not data.exists():
            data.mkdir(parents=True)
            with tarfile.open(ROOT / 'sample_data' / 'sample_1kGP_common_global.tar.gz') as tar:
                for m in tar.getmembers():
                    if m.isfile() and m.name.startswith(('sample_chr1_', 'sample_chr21_', 'sample_chr22_',
                                                         'sample_chrX_')):
                        tar.extract(m, data, filter='data')
        tr = str(ROOT / 'tests' / 'fixtures' / '1kg' / 'transcripts.csv')
        return [(c, str(data / ('sample_chr%s_1kGP_GRCh38_phased*.vcf' % c)), '0.01', str(ROOT / 'igsr_samples.tsv'),
                 tr) for c in ('1', '21', '22', 'X')]

    def simulated(path):
        if not (path / 'transcripts.csv').exists():
            subprocess.run([sys.executable, str(ROOT / 'tests' / 'fixtures' / 'sim' / 'make_data.py'), path.name],
                           check=True, capture_output=True)
        lines = (path / 'transcripts.csv').read_text().split('\n')[1:]
        chroms = list(dict.fromkeys(line.split(',')[0] for line in lines if line))
        return [(c, str(path / ('chr%s*vcf.gz' % c)), '0', str(path / 'samples.tsv'), str(path / 'transcripts.csv'))
                for c in chroms]

    return {'1kg': onekg, 'mix': lambda: simulated(sim / 'mix'), 'ttn': lambda: simulated(sim / 'ttn'),
            'block': lambda: simulated(block)}


def run_prohap(out, chrom, vcf_glob, af, samples, transcripts, fixes=None, script=None):
    """The reference filter, then the fast prohap.py -fix_bugs 1 (PROHAP_FIXES=fixes when given), in out."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    (out / 'tmp').mkdir(exist_ok=True)
    env = dict(os.environ, PYTHONHASHSEED='0', PYTHONNOUSERSITE='1', PROHAP_FALLBACK='error')
    env.pop('PROHAP_FIXES', None)
    if fixes is not None:
        env['PROHAP_FIXES'] = ','.join(fixes)
    runner = [sys.executable] + ([str(ROOT / 'tools' / 'run_fork.py')] if sys.platform == 'darwin' else [])
    filtered = out.parent / ('filtered_chr%s.vcf.gz' % chrom)
    if not filtered.exists():
        subprocess.run(runner + [str(ROOT / 'legacy' / 'src' / 'vcf_filter_fix.py'), '-i', vcf_glob, '-chr', chrom,
                                 '-af', af, '-af_field', 'AF', '-o', str(filtered) + '.part'], check=True,
                       capture_output=True, env=env)
        os.rename(str(filtered) + '.part', filtered)
    db = REFS / ('Homo_sapiens.GRCh38.113.chr_patch_hapl_scaff_chr%s.db' % chrom)
    proc = subprocess.run(runner + [str(script or ROOT / 'src' / 'prohap.py'), '-i', str(filtered), '-db', str(db),
                                    '-transcripts', transcripts, '-cdna', str(REFS / 'total_cdnas_113.fa'),
                                    '-s', samples, '-chr', chrom, '-min_hap_freq', '-1', '-min_hap_count', '0',
                                    '-acc_prefix', 'enshap_' + chrom, '-id_prefix', 'haplo_chr' + chrom,
                                    '-require_start', '1', '-ignore_UTR', '1', '-skip_start_lost', '1',
                                    '-x_par1_to', str(PAR1_TO), '-x_par2_from', str(PAR2_FROM), '-threads', '4',
                                    '-log', 'prohap.log', '-tmp_dir', 'tmp', '-output_csv', 'haplo.tsv',
                                    '-output_fasta', 'haplo.fa'] + ([] if script else ['-fix_bugs', '1']),
                          cwd=out, env=env, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError('prohap.py failed in %s:\n%s' % (out, proc.stderr[-3000:]))
    (out / 'stdout.txt').write_text(proc.stdout)
    return filtered


# ---------------------------------------------------------------------------------------------- inputs

class Transcript:
    def __init__(self, tid, strand, start, end, exons, start_codon, cdna):
        self.id = tid
        self.strand = strand
        self.start = start
        self.end = end
        self.exons = exons                      # [(start, end)] in genome order, inclusive
        self.cdna = cdna
        self.plus = {}                          # genome position -> index in the + strand exonic sequence
        self.exon_of = {}
        for k, (s, e) in enumerate(exons):
            for p in range(s, e + 1):
                self.exon_of[p] = k
                self.plus[p] = len(self.plus)
        self.length = len(self.plus)
        i = self.plus[start_codon]
        self.start_index = i if strand == '+' else self.length - i - 3


def read_transcripts(chrom, transcripts_csv, db_path=None, cdna_path=None):
    import gffutils
    db = gffutils.FeatureDB(str(db_path or REFS / ('Homo_sapiens.GRCh38.113.chr_patch_hapl_scaff_chr%s.db' % chrom)))
    wanted = []
    for line in Path(transcripts_csv).read_text().split('\n')[1:]:
        if line and line.split(',')[0] == chrom:
            wanted.append(line.split(',')[1])
    cdna = {}
    name = None
    with open(cdna_path or REFS / 'total_cdnas_113.fa') as fh:
        keep = set(t.split('.')[0] for t in wanted)
        parts = []
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    cdna[name] = ''.join(parts)
                name = line[1:].split()[0].split('.')[0]
                name = name if name in keep else None
                parts = []
            elif name is not None:
                parts.append(line.strip())
        if name is not None:
            cdna[name] = ''.join(parts)
    result = []
    skipped = []
    for tid in wanted:
        f = db[tid]
        starts = sorted(c.start for c in db.children(f, featuretype='start_codon'))
        if not starts:
            continue                            # require_start
        exons = sorted((e.start, e.end) for e in db.children(f, featuretype='exon'))
        seq = cdna.get(tid.split('.')[0])
        if seq is None:
            skipped.append((tid, 'not in the cDNA FASTA'))
            continue
        t = Transcript(tid, f.strand, f.start, f.end, exons, starts[0], seq)
        if t.length != len(seq):
            skipped.append((tid, 'cDNA is not the spliced exons'))
            continue
        result.append(t)
    return result, skipped


def read_vcf(path, sample_names):
    """records [(pos, id, ref, alt, carriers: set of copy labels)], included sample names in VCF order."""
    records = []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                columns = line.rstrip('\n').split('\t')
                included = [(k, s) for k, s in enumerate(columns) if k >= 9 and s in sample_names]
                continue
            f = line.rstrip('\n').split('\t')
            carriers = set()
            for k, s in included:
                gt = f[k]
                if gt[0] == '1' and gt[1:2] == '|':
                    carriers.add(s + ':1')
                if gt.endswith('|1'):
                    carriers.add(s + ':2')
            records.append((int(f[1]), f[2], f[3], f[4], carriers))
    return records, [s for _, s in included]


def read_metadata(path):
    sex = {}
    lines = Path(path).read_text().split('\n')
    header = lines[0].split('\t')
    for line in lines[1:]:
        if line:
            row = dict(zip(header, line.split('\t')))
            sex[row['Sample name']] = row['Sex']
    return sex


# ---------------------------------------------------------------------------------------------- the rules

def trim(pos, ref, alt):
    """(start, end, ALT) of the bases a record changes; end == start for an insertion before base start."""
    ref = '' if ref == '-' else ref
    alt = '' if alt == '-' else alt
    k = 0
    while k < min(len(ref), len(alt)) and ref[k] == alt[k]:
        k += 1
    ref, alt, pos = ref[k:], alt[k:], pos + k
    m = 0
    while m < min(len(ref), len(alt)) and ref[-1 - m] == alt[-1 - m]:
        m += 1
    return pos, pos + len(ref) - m, alt[:len(alt) - m]


def conflict(changes):
    """Any two changes (start, end, ALT) that overlap, an insertion strictly inside a change, two insertions at
    one point. Pairwise: the oracle is naive on purpose."""
    real = [c for c in changes if not (c[0] == c[1] and c[2] == '')]
    for i in range(len(real)):
        for j in range(i + 1, len(real)):
            (a1, b1, _), (a2, b2, _) = real[i], real[j]
            if a1 < b1 and a2 < b2 and a1 < b2 and a2 < b1:
                return True
            if a1 == b1 and a2 < a1 < b2:
                return True
            if a2 == b2 and a1 < a2 < b1:
                return True
            if a1 == b1 and a2 == b2 and a1 == a2:
                return True
    return False


def expected_protein(t, recs, reference_protein, notes=None):
    """('removed' | 'ref' | protein | ('skip', why)) for one haplotype copy carrying recs [(pos, id, ref, alt)].
    notes, a Counter, gets 'start lost: ...' for a lost start codon: the ATG found again in frame, in another frame,
    or none."""
    changes = [trim(p, r, a) for p, _, r, a in recs]
    if conflict(changes):
        return 'removed'
    plus_seq = t.cdna if t.strand == '+' else rc(t.cdna)
    L = t.length
    kept = []
    for (pos, _, ref, alt), (a, b, alt_t) in zip(recs, changes):
        ref_ = '' if ref == '-' else ref
        span = range(pos, pos + max(len(ref_), 1))
        if any(p not in t.exon_of for p in span) or len({t.exon_of[p] for p in span}) > 1:
            return ('skip', 'record across an exon boundary')
        if a == b and alt_t:
            if (a - 1) not in t.exon_of or a not in t.exon_of or t.exon_of[a - 1] != t.exon_of[a]:
                return ('skip', 'insertion at an exon boundary')
        if plus_seq[t.plus[pos]:t.plus[pos] + len(ref_)] != ref_:
            return ('skip', 'REF not in the cDNA')
        if a == b and not alt_t:
            kept.append((a, b, alt_t))
            continue
        # the liftover check on the trimmed alleles, in the cDNA's orientation
        q = t.plus[a] if a < b else t.plus[a - 1] + 1
        ref_t = plus_seq[q:q + (b - a)]
        if t.strand == '-':
            loc, ref_c, alt_c = L - q - (b - a), rc(ref_t), rc(alt_t)
        else:
            loc, ref_c, alt_c = q, ref_t, alt_t
        if t.cdna[loc:loc + len(alt_c)] == alt_c and t.cdna[loc:loc + len(ref_c)] != ref_c:
            continue
        kept.append((a, b, alt_t))
    if not any(a != b or alt_t for a, b, alt_t in kept):
        return 'ref'
    # the haplotype's exonic sequence (+ strand order), by string edits
    before = collections.defaultdict(str)
    removed = set()
    for a, b, alt_t in sorted(kept, key=lambda c: (c[0], c[1])):
        if a == b:
            before[t.plus[a - 1] + 1] += alt_t
        else:
            for p in range(a, b):
                removed.add(t.plus[p])
            before[t.plus[a]] += alt_t
    plus_hap = ''.join(before[q] + ('' if q in removed else plus_seq[q]) for q in range(L)) + before[L]
    hap = plus_hap if t.strand == '+' else rc(plus_hap)
    # the cDNA intervals of the changes, to place the start codon
    intervals = []
    for a, b, alt_t in kept:
        if a == b:
            q = t.plus[a - 1] + 1
            lo, hi = (q, q) if t.strand == '+' else (L - q, L - q)
        else:
            lo, hi = (t.plus[a], t.plus[b - 1] + 1) if t.strand == '+' else (L - t.plus[b - 1] - 1, L - t.plus[a])
        intervals.append((lo, hi, len(alt_t) - (b - a)))
    s = t.start_index

    def shift(x):
        return sum(d for lo, hi, d in intervals if hi <= x)
    touching = [(lo, hi) for lo, hi, d in intervals
                if (lo < hi and lo < s + 3 and s < hi) or (lo == hi and s < lo < s + 3)]
    if not touching:
        return translate(hap[s + shift(s):])
    p0 = min([s] + [lo for lo, hi in touching])
    q = hap.find('ATG', p0 + shift(p0))
    if notes is not None:
        notes['start lost: ' + ('no ATG' if q == -1 else 'ATG in frame' if (q - s - shift(s)) % 3 == 0
                                else 'ATG in another frame')] += 1
    if q == -1:
        return 'ref'
    return translate(hap[q:])


# ---------------------------------------------------------------------------------------------- comparison

def fast_outcome(out):
    sequences = {}
    with open(out / 'haplo.fa') as fh:
        lines = fh.read().split('\n')
    for header, seq in zip(lines[0::2], lines[1::2]):
        for h in header.split('matching_proteins:')[1].split()[0].split(';'):
            sequences[h] = seq
    result = {}
    rows = {}
    with open(out / 'haplo.tsv') as fh:
        columns = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            row = dict(zip(columns, line.rstrip('\n').split('\t')))
            for copy in row['samples'].split(';'):
                result[(row['TranscriptID'], copy)] = sequences.get(row['HaplotypeID'])
                rows[(row['TranscriptID'], copy)] = row
    with open(out / 'prohap.log') as fh:
        fh.readline()
        for line in fh:
            transcript, _, removed = line.rstrip('\n').split('\t')
            for copy in (removed.split(';') if removed else []):
                result[(transcript, copy)] = 'removed'
    return result, rows


def compare(out, filtered, chrom, samples_tsv, transcripts_csv, examples=5, db_path=None, cdna_path=None):
    """Counts of the comparison for one chromosome, and examples of each kind of difference."""
    out = Path(out)
    sex = read_metadata(samples_tsv)
    records, included = read_vcf(filtered, set(sex))
    transcripts, skipped = read_transcripts(chrom, transcripts_csv, db_path, cdna_path)
    got, rows = fast_outcome(out)
    counts = collections.Counter()
    shown = collections.defaultdict(list)
    by_pos = sorted(range(len(records)), key=lambda i: records[i][0])
    for t in transcripts:
        recs = [records[i] for i in by_pos
                if any(s < records[i][0] + len(records[i][2]) and e >= records[i][0] for s, e in t.exons)]
        if not recs:
            continue
        autosomal = chrom != 'X' or (t.start < PAR1_TO and t.end <= PAR1_TO) or \
            (t.start >= PAR2_FROM and t.end > PAR2_FROM)
        copies = [s + h for s in included for h in (':1', ':2') if autosomal or h == ':1' or sex[s] != 'male']
        reference_protein = translate(t.cdna[t.start_index:])
        id_count = collections.Counter('NA' if r[1] in NA_TOKENS else r[1] for r in recs if r[1] != '.')
        carried = collections.defaultdict(list)
        for k, r in enumerate(recs):
            for c in r[4]:
                carried[c].append(k)
        by_set = collections.defaultdict(list)
        for c in copies:
            by_set[tuple(carried.get(c, ()))].append(c)
        for key, members in by_set.items():
            hap = [recs[k][:4] for k in key]
            notes = collections.Counter()
            want = expected_protein(t, hap, reference_protein, notes) if hap else 'ref'
            for note in notes:
                counts[note + ' (haplotypes)'] += 1
                counts[note + ' (copies)'] += len(members)
            if want not in ('removed', 'ref') and not isinstance(want, tuple) and want == reference_protein:
                want = 'ref'
            for c in members:
                have = got.get((t.id, c), 'ref')
                if isinstance(want, tuple):
                    counts['left out: ' + want[1]] += 1
                    continue
                counts['copies compared'] += 1
                if have == want:
                    counts['same'] += 1
                    continue
                row = rows.get((t.id, c))
                kind = 'DIFFERENT'
                if row is not None and any(id_count['NA' if i in NA_TOKENS else i] > 1
                                           for i in row['VCF_IDs'].split(';')):
                    kind = 'shared_id'
                counts[kind] += 1
                if len(shown[kind]) < examples:
                    shown[kind].append('%s %s %s: fast %s, oracle %s' % (t.id, t.strand, c, short(have), short(want))
                                       + ' records ' + ' '.join('%d:%s>%s' % (p, r, a) for p, _, r, a in hap))
    for tid, why in skipped:
        counts['transcripts left out: ' + why] += 1
    return counts, shown


def short(x):
    return x if len(x) <= 40 else x[:18] + '...' + x[-18:] + ' (%d aa)' % len(x)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('datasets', nargs='+', help='1kg, mix, ttn, block, or all')
    parser.add_argument('--work', help='work directory (default: a temporary one)')
    parser.add_argument('--keep', action='store_true')
    parser.add_argument('--chromosomes', help='only these chromosomes, comma-separated')
    args = parser.parse_args(argv)
    work = Path(args.work or tempfile.mkdtemp(prefix='fix_oracle.')).resolve()
    work.mkdir(parents=True, exist_ok=True)
    known = datasets(work)
    names = list(known) if args.datasets == ['all'] else args.datasets
    total = collections.Counter()
    for name in names:
        for chrom, vcf_glob, af, samples, transcripts in known[name]():
            if args.chromosomes and chrom not in args.chromosomes.split(','):
                continue
            out = work / name / ('chr' + chrom) / 'fast'
            filtered = run_prohap(out, chrom, vcf_glob, af, samples, transcripts)
            counts, shown = compare(out, filtered, chrom, samples, transcripts)
            total.update(counts)
            print('%s chr%s: %s' % (name, chrom, ', '.join('%s %d' % kv for kv in sorted(counts.items()))))
            for kind, lines in shown.items():
                for line in lines:
                    print('   %s: %s' % (kind, line))
    print('TOTAL: ' + ', '.join('%s %d' % kv for kv in sorted(total.items())))
    if not args.keep and not args.work:
        shutil.rmtree(work)
    return 1 if total['DIFFERENT'] else 0


if __name__ == '__main__':
    sys.exit(main())
