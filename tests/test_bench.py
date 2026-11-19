"""bench/simulate.py, bench/simulate_fasta.py and bench/bench.py on a tiny synthetic references directory, and the
consistency of the simulated fixtures (tests/fixtures/sim) with tests/legacy_vs_fast.sh.

Runs in the ProHap environment, without pre-built references (about 30 s):

    pytest tests/test_bench.py
"""
import gzip
import json
import random
import re
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip('gffutils')
pytest.importorskip('numpy')

ROOT = Path(__file__).resolve().parent.parent
SIMULATE = ROOT / 'bench' / 'simulate.py'
SIMULATE_FASTA = ROOT / 'bench' / 'simulate_fasta.py'
BENCH = ROOT / 'bench' / 'bench.py'
SIM = ROOT / 'tests' / 'fixtures' / 'sim'
RELEASE = '113'
PAR1_TO, PAR2_FROM = 2781479, 155701383
COMPLEMENT = str.maketrans('ACGT', 'TGCA')

# (chromosome, transcript, strand, exons, 5' UTR length, codons from ATG to the stop, start codon annotated)
TRANSCRIPTS = [
    ('22', 'ENST00000000001', '+', [(10000, 10599), (11000, 11999), (13000, 13799)], 50, 600, True),
    ('22', 'ENST00000000002', '-', [(20000, 20999), (22000, 22599)], 40, 450, True),
    ('22', 'ENST00000000003', '+', [(30000, 30299)], 10, 60, False),        # no start codon: never in a block
    ('X', 'ENST00000000004', '+', [(100000, 100899)], 30, 250, True),        # PAR1
    ('X', 'ENST00000000005', '-', [(3000000, 3001199)], 60, 300, True),
    ('X', 'ENST00000000006', '+', [(3100000, 3100599)], 20, 150, True),
]


def genomic_positions(strand, exons):
    """Genomic position of every cDNA index."""
    positions = [p for s, e in exons for p in range(s, e + 1)]
    return positions if strand == '+' else positions[::-1]


def build_transcript(rng, strand, exons, utr5, codons):
    length = sum(e - s + 1 for s, e in exons)
    sense = [c for c in (a + b + c for a in 'ACGT' for b in 'ACGT' for c in 'ACGT') if c not in ('TAA', 'TAG', 'TGA')]
    cds = 'ATG' + ''.join(rng.choice(sense) for _ in range(codons - 2)) + 'TAA'
    cdna = ''.join(rng.choice('ACGT') for _ in range(utr5)) + cds
    cdna += ''.join(rng.choice('ACGT') for _ in range(length - len(cdna)))
    pos = genomic_positions(strand, exons)
    codon = sorted(pos[utr5:utr5 + 3])
    return cdna, (codon[0], codon[-1])


@pytest.fixture(scope='module')
def refs(tmp_path_factory):
    """A references directory with the files simulate.py, simulate_fasta.py and prohap.py read."""
    import gffutils
    rng = random.Random(1)
    out = tmp_path_factory.mktemp('refs') / ('ensembl.v' + RELEASE)
    out.mkdir()
    (out / 'release.txt').write_text(RELEASE + '\n')
    gtf = {}
    cdna_fa = []
    pep_fa = []
    for chrom, tid, strand, exons, utr5, codons, has_start in TRANSCRIPTS:
        cdna, start = build_transcript(rng, strand, exons, utr5, codons)
        gid = tid.replace('ENST', 'ENSG')
        attrs = 'gene_id "%s"; transcript_id "%s"; transcript_biotype "protein_coding";' % (gid, tid)
        lines = gtf.setdefault(chrom, [])
        lines.append('\t'.join([chrom, 'test', 'transcript', str(exons[0][0]), str(exons[-1][1]), '.', strand, '.',
                                attrs]))
        lines += ['\t'.join([chrom, 'test', 'exon', str(s), str(e), '.', strand, '.', attrs]) for s, e in exons]
        if has_start:
            lines.append('\t'.join([chrom, 'test', 'start_codon', str(start[0]), str(start[1]), '.', strand, '0',
                                    attrs]))
        cdna_fa.append('>%s.1 cdna chromosome:GRCh38:%s:%d:%d:1\n%s\n' % (tid, chrom, exons[0][0], exons[-1][1], cdna))
        pep_fa.append('>generic_ensref|%s|pep chromosome:GRCh38:%s:%d:%d:1 matching_proteins:%s\n%s\n'
                      % (tid, chrom, exons[0][0], exons[-1][1], tid, 'M' + 'ACDEFGHIKL' * 20))
    for chrom, lines in gtf.items():
        path = out / ('chr%s.gtf' % chrom)
        path.write_text('\n'.join(lines) + '\n')
        gffutils.create_db(str(path), dbfn=str(out / ('Homo_sapiens.GRCh38.%s.chr_patch_hapl_scaff_chr%s.db'
                                                      % (RELEASE, chrom))),
                           force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True,
                           disable_infer_genes=True, disable_infer_transcripts=True)
    (out / ('total_cdnas_%s.fa' % RELEASE)).write_text(''.join(cdna_fa))
    (out / ('ensembl_reference_proteinDB_%s_tagged.fa' % RELEASE)).write_text(''.join(pep_fa))
    (out / ('transcripts_reference_%s.csv' % RELEASE)).write_text(
        'chromosome,transcriptID\n' + ''.join('%s,%s\n' % (t[0], t[1]) for t in TRANSCRIPTS) + '7,ENST00000000099\n')
    with gzip.open(out / 'crap_tagged.fa.gz', 'wt') as fh:
        fh.write('>generic_sp|ALBU_BOVIN| matching_proteins:ALBU_BOVIN\nMKWVTFISLLLLFSSAYSRGVFRR\n')
    return out


def genome(refs_dir):
    """{(chromosome, position): + strand base} from the tiny references."""
    cdnas = {}
    for block in (refs_dir / ('total_cdnas_%s.fa' % RELEASE)).read_text().split('>')[1:]:
        header, seq = block.split('\n', 1)
        cdnas[header.split('.')[0]] = seq.strip()
    bases = {}
    for chrom, tid, strand, exons, *_ in TRANSCRIPTS:
        for p, b in zip(genomic_positions(strand, exons), cdnas[tid]):
            bases[(chrom, p)] = b if strand == '+' else b.translate(COMPLEMENT)
    return bases


def simulate(refs_dir, out, *args, check=True):
    return subprocess.run([sys.executable, str(SIMULATE), '--refs', str(refs_dir), '--out', str(out)] + list(args),
                          capture_output=True, text=True, check=check)


def records(path):
    with gzip.open(path, 'rt') as fh:
        lines = fh.read().split('\n')
    assert lines[-1] == ''
    header = [line for line in lines if line.startswith('#')]
    return header, [line.split('\t') for line in lines[:-1] if not line.startswith('#')]


@pytest.fixture(scope='module')
def dataset(refs, tmp_path_factory):
    out = tmp_path_factory.mktemp('sim')
    simulate(refs, out, '--chrom', '22', '--transcripts', '0', '--samples', '2000', '--ids', 'dot', '--unlisted', '3')
    simulate(refs, out, '--chrom', 'X', '--transcripts', '0', '--samples', '2000', '--ids', 'rs', '--unlisted', '3')
    return out


def test_same_arguments_same_bytes(refs, dataset, tmp_path):
    simulate(refs, tmp_path, '--chrom', '22', '--transcripts', '0', '--samples', '2000', '--ids', 'dot',
             '--unlisted', '3')
    for name in ('chr22.vcf.gz', 'samples.tsv'):
        assert (tmp_path / name).read_bytes() == (dataset / name).read_bytes(), name


def test_block_skips_transcripts_without_start_codon(dataset):
    assert (dataset / 'transcripts.csv').read_text() == (
        'chromosome,transcriptID\n22,ENST00000000001\n22,ENST00000000002\n'
        'X,ENST00000000004\nX,ENST00000000005\nX,ENST00000000006\n')
    summary = json.loads((dataset / 'chr22.json').read_text())
    assert summary['block'] == ['ENST00000000001', 'ENST00000000002']
    assert summary['exonic_bp'] == 4000


def test_ref_matches_the_reference(refs, dataset):
    bases = genome(refs)
    _, rows = records(dataset / 'chr22.vcf.gz')
    assert len(rows) > 50
    positions = [int(r[1]) for r in rows]
    assert positions == sorted(positions)
    kinds = set()
    for r in rows:
        pos, ref, alts = int(r[1]), r[3], r[4].split(',')
        assert re.fullmatch('[ACGT]+', ref) and all(re.fullmatch('[ACGT]+', a) for a in alts)
        # every REF base on an exon is the reference base; bases past the exon are not checked by anything
        for k, b in enumerate(ref):
            if ('22', pos + k) in bases:
                assert bases[('22', pos + k)] == b, (pos, k)
        assert ('22', pos) in bases
        if len(ref) > 1:
            kinds.add('deletion')
            assert alts == [ref[0]]
        elif len(alts[0]) > 1:
            kinds.add('insertion')
            assert alts[0][0] == ref
        else:
            kinds.add('multiallelic' if len(alts) > 1 else 'snv')
            assert ref not in alts and len(set(alts)) == len(alts)
    assert kinds == {'snv', 'deletion', 'insertion', 'multiallelic'}


def test_genotypes_and_info(dataset):
    header, rows = records(dataset / 'chr22.vcf.gz')
    names = header[-1].split('\t')[9:]
    assert header[-1].startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
    assert len(names) == 2000 and names[0] == 'S000001'
    for r in rows:
        assert r[2] == '.' and r[5] == '.' and r[6] == 'PASS' and r[8] == 'GT'
        gts = r[9:]
        assert len(gts) == 2000
        n_alts = len(r[4].split(','))
        alleles = [a for gt in gts for a in gt.split('|')]
        assert all(len(gt) == 3 and gt[1] == '|' for gt in gts)
        ac = [alleles.count(str(k)) for k in range(1, n_alts + 1)]
        info = dict(kv.split('=') for kv in r[7].split(';'))
        assert info['AC'] == ','.join(map(str, ac)) and info['AN'] == '4000'
        assert info['AF'] == ','.join('%.6g' % (a / 4000) for a in ac)
        assert all(a >= 1 for a in ac)


def test_chrx_males(refs, dataset, tmp_path):
    samples = (dataset / 'samples.tsv').read_text().split('\n')[1:-1]
    male = {line.split('\t')[0] for line in samples if line.split('\t')[1] == 'male'}
    header, haploid = records(dataset / 'chrX.vcf.gz')
    names = header[-1].split('\t')[9:]
    simulate(refs, tmp_path, '--chrom', 'X', '--transcripts', '0', '--samples', '2000', '--ids', 'rs',
             '--unlisted', '3', '--males', 'mirrored')
    _, mirrored = records(tmp_path / 'chrX.vcf.gz')
    assert [r[:7] for r in haploid] == [r[:7] for r in mirrored]
    seen_par = seen_non_par = False
    for h, m in zip(haploid, mirrored):
        non_par = PAR1_TO < int(h[1]) < PAR2_FROM
        seen_par |= not non_par
        seen_non_par |= non_par
        haploid_tokens = 0
        for name, gh, gm in zip(names, h[9:], m[9:]):
            if non_par and len(gh) == 1:
                haploid_tokens += 1
                assert gm == gh + '|' + gh
            else:
                assert len(gh) == 3 and gh == gm
            # the metadata leaves out 3 of the VCF's samples (--unlisted)
            assert name not in male or len(gh) == (1 if non_par else 3)
        if non_par:
            assert len(male) - 3 <= haploid_tokens <= len(male) + 3
        for row, an in ((h, 4000 - haploid_tokens), (m, 4000)):
            info = dict(kv.split('=') for kv in row[7].split(';'))
            ac = [sum(g.split('|').count(str(k)) if '|' in g else (g == str(k)) for g in row[9:])
                  for k in range(1, len(row[4].split(',')) + 1)]
            assert info['AN'] == str(an) and info['AC'] == ','.join(map(str, ac))
    assert seen_par and seen_non_par


def test_ids(refs, dataset, tmp_path):
    _, rs = records(dataset / 'chrX.vcf.gz')
    assert all(re.fullmatch(r'rs[0-9]+', r[2]) for r in rs) and len({r[2] for r in rs}) == len(rs)
    simulate(refs, tmp_path, '--chrom', '22', '--transcripts', '0', '--samples', '2000', '--ids', 'ambiguous',
             '--unlisted', '3')
    _, ambiguous = records(tmp_path / 'chr22.vcf.gz')
    _, dot = records(dataset / 'chr22.vcf.gz')
    assert [r[:2] + r[3:] for r in ambiguous] == [r[:2] + r[3:] for r in dot]     # only the IDs differ
    tokens = [r[2] for r in ambiguous if r[2] != '.']
    assert 0.15 < len(tokens) / len(ambiguous) < 0.45
    assert all(re.fullmatch(r'[0-9]e[0-9]|[0-9]+|00[0-9]+|NA|nan|True|False|[0-9]\.[0-9]0', t) for t in tokens)


def test_chunks(refs, dataset, tmp_path):
    simulate(refs, tmp_path, '--chrom', '22', '--transcripts', '0', '--samples', '2000', '--unlisted', '3',
             '--chunks', '2')
    assert sorted(p.name for p in tmp_path.glob('chr22*.vcf.gz')) == ['chr22.1.vcf.gz', 'chr22.2.vcf.gz']
    header1, part1 = records(tmp_path / 'chr22.1.vcf.gz')
    header2, part2 = records(tmp_path / 'chr22.2.vcf.gz')
    whole_header, whole = records(dataset / 'chr22.vcf.gz')
    assert header1 == header2 and header1[-1] == whole_header[-1]
    assert part1 and part2 and part1 + part2 == whole
    simulate(refs, tmp_path, '--chrom', '22', '--transcripts', '0', '--samples', '2000', '--unlisted', '3')
    assert sorted(p.name for p in tmp_path.glob('chr22*.vcf.gz')) == ['chr22.vcf.gz']


def test_samples(refs, dataset, tmp_path):
    lines = (dataset / 'samples.tsv').read_text().split('\n')
    assert lines[0] == 'Sample name\tSex\tPopulation code\tPopulation name\tSuperpopulation code\tSuperpopulation name'
    rows = [line.split('\t') for line in lines[1:-1]]
    header, _ = records(dataset / 'chr22.vcf.gz')
    vcf = header[-1].split('\t')[9:]
    listed = [r[0] for r in rows]
    assert len(listed) == 2000 and len(set(vcf) - set(listed)) == 3
    assert [n for n in listed if n not in vcf] == ['M000001', 'M000002', 'M000003']
    assert {r[1] for r in rows} == {'male', 'female'} and len({r[2] for r in rows}) == 26
    result = simulate(refs, tmp_path / 'other', '--chrom', '22', '--transcripts', '0', '--samples', '2000',
                      check=False)
    assert result.returncode == 0
    result = simulate(refs, dataset, '--chrom', '22', '--transcripts', '0', '--samples', '1999',
                      check=False)
    assert result.returncode != 0 and 'exists with other samples' in result.stderr


def run_reference(refs, dataset, chrom, work):
    """The legacy filter and prohap.py, as the workflow runs them; returns prohap.py's stdout."""
    runner = [sys.executable] + ([str(ROOT / 'tools' / 'run_fork.py')] if sys.platform == 'darwin' else [])
    env = {'PYTHONHASHSEED': '0', 'PYTHONNOUSERSITE': '1', 'PATH': '/usr/bin:/bin'}
    work.mkdir()
    (work / 'tmp').mkdir()
    subprocess.run(runner + [str(ROOT / 'legacy' / 'src' / 'vcf_filter_fix.py'), '-i',
                             str(dataset / ('chr%s.*vcf.gz' % chrom)), '-chr', chrom, '-af', '0', '-o',
                             'filtered.vcf.gz'], cwd=work, env=env, check=True, capture_output=True)
    return subprocess.run(runner + [str(ROOT / 'legacy' / 'src' / 'prohap.py'), '-i', 'filtered.vcf.gz', '-db',
                                    str(refs / ('Homo_sapiens.GRCh38.%s.chr_patch_hapl_scaff_chr%s.db'
                                                % (RELEASE, chrom))),
                                    '-transcripts', str(dataset / 'transcripts.csv'),
                                    '-cdna', str(refs / ('total_cdnas_%s.fa' % RELEASE)),
                                    '-s', str(dataset / 'samples.tsv'), '-chr', chrom, '-min_hap_freq', '0',
                                    '-min_hap_count', '0', '-threads', '1', '-log', 'removed.log', '-tmp_dir', 'tmp',
                                    '-output_csv', 'haplo.tsv', '-output_fasta', 'haplo.fa'],
                          cwd=work, env=env, check=True, capture_output=True, text=True).stdout


@pytest.mark.parametrize('chrom', ['22', 'X'])
def test_reference_code_reads_the_data(refs, dataset, tmp_path, chrom):
    out = run_reference(refs, dataset, chrom, tmp_path / 'run')
    assert 'Ref allele not matching' not in out and 'Incorrect formatting' not in out
    table = (tmp_path / 'run' / 'haplo.tsv').read_text().split('\n')
    assert len(table) > 10
    assert (tmp_path / 'run' / 'haplo.fa').read_text().count('>') > 5


def test_bench_harness(refs, dataset, tmp_path):
    tsv = tmp_path / 'results.tsv'
    result = subprocess.run([sys.executable, str(BENCH), '--refs', str(refs), '--data', str(dataset), '--arms',
                             'legacy,fast', '--tsv', str(tsv), '--work', str(tmp_path / 'work')],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'dataset chr22, legacy vs fast: outputs IDENTICAL' in result.stdout.replace(dataset.name, 'dataset')
    rows = [line.split('\t') for line in tsv.read_text().split('\n')[1:-1]]
    assert [(r[3], r[7], r[8]) for r in rows] == [(c, a, s) for c in ('22', 'X') for a in ('legacy', 'fast')
                                                  for s in ('filter', 'prohap')]
    prohap = [r for r in rows if r[8] == 'prohap']
    assert all(float(r[11]) > 0 and float(r[13]) > 0 for r in prohap)
    expected = {'legacy': ['setup', 'sweep', 'haplotypes', 'cdna', 'proteins', 'write'],
                'fast': ['setup', 'cdna', 'sweep', 'write']}
    assert all([p.split('=')[0] for p in r[14].split(';')] == expected[r[7]] for r in prohap)


def test_simulate_fasta(refs, tmp_path):
    def run(out):
        subprocess.run([sys.executable, str(SIMULATE_FASTA), '--refs', str(refs), '--records', '300', '--out',
                        str(out)], check=True, capture_output=True)
        return gzip.open(out, 'rt').read()
    text = run(tmp_path / 'a.fa.gz')
    assert text == run(tmp_path / 'b.fa.gz')
    headers = [line for line in text.split('\n') if line.startswith('>')]
    assert len(headers) == 300
    assert headers[0].startswith('>generic_ensref|ENST00000000001|position_within_protein:0 pep ')
    assert headers[6] == '>generic_sp|ALBU_BOVIN| matching_proteins:ALBU_BOVIN'
    assert all(re.fullmatch(r'>generic_enshap\|enshap_(22|X)_[0-9a-f]+(_3UTR_1)?\|position_within_protein:[0-9]+ '
                            r'matching_proteins:[^ ]+ start:0 reading_frame:0(;0)*', h) for h in headers[7:])
    sequences = text.split('\n')[1::2]
    assert len(set(sequences)) < len(sequences)     # duplicates for merge_duplicate_seq.py


def test_fixture_definitions():
    """Every simulated config names a dataset of make_data.py with the same chromosomes, and every tier of
    legacy_vs_fast.sh names an existing fixture."""
    sys.path.insert(0, str(SIM))
    import make_data
    configs = sorted(SIM.glob('*.yaml'))
    assert configs
    for config in configs:
        dataset = make_data.DATASETS[config.stem.split('-')[0]]
        chroms = re.search(r'^chromosomes: \[(.*)\]$', config.read_text(), re.M).group(1)
        assert sorted(c.strip(' "') for c in chroms.split(',')) == \
            sorted(args[args.index('--chrom') + 1] for args in dataset['chromosomes']), config.name
    script = (ROOT / 'tests' / 'legacy_vs_fast.sh').read_text()
    named = set()
    for tier in ('SMALL', 'MEDIUM', 'LARGE'):
        named.update(re.search(r'^%s="(.*)"$' % tier, script, re.M).group(1).split())
    for fixture in named:
        assert fixture == 'dupseq' or (ROOT / 'tests' / 'fixtures' / (fixture + '.yaml')).is_file(), fixture
    assert {'sim/' + c.stem for c in configs} | {'1kg/threshold0', '1kg/template', 'dupseq'} == named
    assert set(make_data.DATASETS) == {c.stem.split('-')[0] for c in configs} | {'dupseq'}
