"""verify's preparation: config rewriting, subsetting and the preflight checks (verify/verify.py, verify/vcfcut.py,
verify/preflight.py).

    pytest tests/test_verify_prepare.py

The dry run needs snakemake and gffutils, from the ProHap environment; it is skipped without them.
"""
import argparse
import gzip
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'verify'))
sys.path.insert(0, str(ROOT / 'tests'))
import preflight  # noqa: E402
import vcfcut  # noqa: E402
import verify  # noqa: E402

BIN = Path(sys.executable).parent

HEADER = ('##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description="All filters passed">\n##contig=<ID=22>\n'
          '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n'
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')


def vcf_text(samples, records, header=HEADER):
    """records: (pos, id, genotypes) with one genotype per sample, or a full line ending in a newline."""
    out = header + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples) + '\n'
    for r in records:
        if isinstance(r, str):
            out += r
        else:
            pos, rid, gts = r
            out += '\t'.join(['22', str(pos), rid, 'A', 'G', '.', 'PASS', 'AF=0.0123456789', 'GT'] + gts) + '\n'
    return out


def write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.name.endswith('.gz'):
        with gzip.open(path, 'wt') as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


def read(path):
    path = str(path)
    with (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)) as fh:
        return fh.read()


def args(**kw):
    base = dict(cores=1, chromosomes=None, samples=None, sample_list=None, seed=0, transcripts=None)
    base.update(kw)
    return argparse.Namespace(**base)


def config(**kw):
    cfg = {
        'implementation': 'fast', 'ensembl_release': 113, 'refs_dir': 'refs', 'chromosomes': ['21', '22'],
        'contaminants_fasta': 'crap.fasta', 'only_MANE_select': True, 'final_fasta_file': 'results/db.fa',
        'use_ProHap': True, 'phased_local_path': 'vcfs/', 'phased_vcf_file_name': 'chr{chr}*.vcf',
        'sample_metadata_file': 'samples.tsv', 'haplo_fasta_file': 'results/haplo.fa', 'haplo_cdna_file': '',
        'haplo_table_file': 'results/haplo.tsv.gz', 'use_ProVar': False, 'variant_vcf': {},
        'var_fasta_file': '', 'var_cdna_file': '', 'var_table_file': '', 'haplo_added_table': '',
        'haplo_added_fasta': '',
    }
    cfg.update(kw)
    return cfg


# ------------------------------------ config rewriting ------------------------------------

def test_arm_config_resolves_inputs_and_moves_outputs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = config(final_fasta_file='/abs/out/db.fa', haplo_table_file='../elsewhere/haplo.tsv.gz',
                 haplo_added_fasta='/data/added.fa', custom_transcript_list='lists/tr.csv',
                 variant_vcf={'gnomad': {'file': 'vcf/gnomad.vcf', 'fasta_accession_prefix': 'var', 'min_af': 0}},
                 chromosomes=[21, 22])
    run = verify.Run(cfg, args(), str(tmp_path))
    legacy, fast = run.arm_config('legacy'), run.arm_config('fast')
    assert legacy['implementation'] == 'legacy' and fast['implementation'] == 'fast'
    assert {k: v for k, v in legacy.items() if k != 'implementation'} == \
        {k: v for k, v in fast.items() if k != 'implementation'}
    up = os.path.join('..', '..', '..')
    assert fast['refs_dir'] == str(tmp_path / 'refs')
    assert fast['contaminants_fasta'] == os.path.join(up, 'crap.fasta')
    assert fast['sample_metadata_file'] == os.path.join(up, 'samples.tsv')
    assert fast['phased_local_path'] == os.path.join(up, 'vcfs/')
    assert fast['phased_vcf_file_name'] == 'chr{chr}*.vcf'
    assert fast['custom_transcript_list'] == os.path.join(up, 'lists/tr.csv')
    assert fast['variant_vcf']['gnomad'] == {'file': os.path.join(up, 'vcf/gnomad.vcf'),
                                            'fasta_accession_prefix': 'var', 'min_af': 0}
    assert fast['haplo_added_fasta'] == '/data/added.fa'
    assert fast['haplo_added_table'] == ''
    # relative outputs stay; absolute ones and those leaving the work directory move under outputs/
    assert fast['haplo_fasta_file'] == 'results/haplo.fa'
    assert fast['final_fasta_file'] == 'outputs/abs/out/db.fa'
    assert fast['haplo_table_file'] == os.path.join('outputs', str(tmp_path.parent / 'elsewhere/haplo.tsv.gz')[1:])
    assert fast['haplo_cdna_file'] == ''
    assert fast['chromosomes'] == [21, 22]
    # an arm's relative inputs resolve to the same files as the user's
    arm = tmp_path / 'verify_runs' / 'x' / 'fast'
    assert os.path.normpath(arm / fast['contaminants_fasta']) == str(tmp_path / 'crap.fasta')


def test_empty_phased_local_path_means_data_vcf_phased(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    run = verify.Run(config(phased_local_path=''), args(), str(tmp_path))
    assert run.pattern == 'data/vcf/phased/chr{chr}*.vcf'
    assert run.arm_config('fast')['phased_local_path'] == os.path.join('..', '..', '..', 'data/vcf/phased/')


def test_refs_dir_defaults_like_the_snakefile(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('PROHAP_REFS', 'my_refs')
    assert verify.refs_root(config(refs_dir='')) == str(tmp_path / 'my_refs')
    monkeypatch.delenv('PROHAP_REFS')
    assert verify.refs_root(config(refs_dir='')) == '/opt/prohap/refs'


@pytest.mark.parametrize('custom,mane,expected', [
    (None, True, ('file', 'R/ensembl.v113/transcripts_reference_113.csv')),
    (None, False, ('rule', 'data/transcripts_reference_113.csv')),
    ('data/transcripts_reference_113.csv', True, ('file', 'R/ensembl.v113/transcripts_reference_113.csv')),
    ('./data/transcripts_reference_113.csv', False, ('rule', 'data/transcripts_reference_113.csv')),
    ('', True, ('rule', 'data/included_transcripts.csv')),
    ('mine.csv', True, ('file', 'mine.csv')),
])
def test_transcript_source_mirrors_the_snakefile(custom, mane, expected):
    cfg = config(only_MANE_select=mane)
    if custom is not None:
        cfg['custom_transcript_list'] = custom
    assert verify.transcript_source(cfg, 'R') == expected


def test_workdir_reference_list_is_not_resolved(tmp_path, monkeypatch):
    """The Snakefile treats the upstream value data/transcripts_reference_<rel>.csv specially; keep it verbatim."""
    monkeypatch.chdir(tmp_path)
    cfg = config(custom_transcript_list='data/transcripts_reference_113.csv')
    assert verify.Run(cfg, args(), str(tmp_path)).arm_config('fast')['custom_transcript_list'] == \
        'data/transcripts_reference_113.csv'


def test_chromosome_subset_keeps_config_order(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    run = verify.Run(config(chromosomes=[1, 2, 'X', 22]), args(chromosomes='22,X'), str(tmp_path))
    assert run.chromosomes == ['X', 22]
    with pytest.raises(verify.VerifyError, match='not in the config'):
        verify.Run(config(), args(chromosomes='3'), str(tmp_path))


# ------------------------------------ subsetting ------------------------------------

SAMPLES = ['S%02d' % i for i in range(40)]


def one_vcf(tmp_path):
    write(tmp_path / 'vcfs' / 'chr22.vcf', vcf_text(SAMPLES, [(100, '.', ['0|1'] * 40)]))


def test_samples_are_seeded_and_in_header_order(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    one_vcf(tmp_path)
    pick = lambda seed: verify.select_samples(verify.Run(config(), args(samples=10, seed=seed), str(tmp_path)))
    a = pick(0)
    assert len(a) == 10 and a == sorted(a, key=SAMPLES.index)
    assert pick(0) == a
    assert pick(1) != a
    assert verify.select_samples(verify.Run(config(), args(samples=40), str(tmp_path))) is None
    with pytest.raises(verify.VerifyError, match='only 40 samples'):
        pick_41 = verify.Run(config(), args(samples=41), str(tmp_path))
        verify.select_samples(pick_41)


def test_sample_list(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    one_vcf(tmp_path)
    (tmp_path / 'ids.txt').write_text('S30\n\nS02\n# comment\nS10\n')
    run = verify.Run(config(), args(sample_list='ids.txt'), str(tmp_path))
    assert verify.select_samples(run) == ['S02', 'S10', 'S30']
    (tmp_path / 'bad.txt').write_text('S02\nNOPE\n')
    with pytest.raises(verify.VerifyError, match='1 samples are not in'):
        verify.select_samples(verify.Run(config(), args(sample_list='bad.txt'), str(tmp_path)))


def test_subset_vcfs_keeps_every_record_and_chunk_names(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    gts = lambda k: ['%d|%d' % ((i + k) % 2, (i * k) % 2) for i in range(40)]
    records = [(100, '.', gts(1)), (200, 'rs1', gts(2)), (200, 'rs2', ['1'] * 20 + ['0|1'] * 20),
               (300, '.', ['.|0'] * 40)]
    write(tmp_path / 'vcfs' / 'chr22_a.vcf.gz', vcf_text(SAMPLES, records[:2]))
    write(tmp_path / 'vcfs' / 'chr22_b.vcf.gz', vcf_text(SAMPLES, records[2:]))
    write(tmp_path / 'vcfs' / 'chr21.vcf', vcf_text(SAMPLES, records))
    run = verify.Run(config(phased_vcf_file_name='chr{chr}*.vcf*'), args(samples=7), str(tmp_path))
    os.makedirs(run.input_dir)
    samples = verify.select_samples(run)
    verify.subset_vcfs(run, samples)
    assert (Path(run.input_dir) / 'samples.txt').read_text() == ''.join(s + '\n' for s in samples)
    cols = [9 + SAMPLES.index(s) for s in samples]
    for chrom, names in (('22', ['chr22_a.vcf.gz', 'chr22_b.vcf.gz']), ('21', ['chr21.vcf'])):
        out = Path(run.input_dir) / 'vcf' / ('chr' + chrom)
        assert sorted(os.listdir(out)) == names
        for name in names:
            original, cut = read(tmp_path / 'vcfs' / name), read(out / name)
            if name.endswith('.gz'):
                assert (out / name).read_bytes()[:2] == b'\x1f\x8b'
            # only the other samples' columns are gone
            expected = [l if l.startswith('##') else '\t'.join([f for i, f in enumerate(l.split('\t'))
                                                                 if i < 9 or i in cols])
                        for l in original.split('\n')]
            assert cut.split('\n') == expected
    fast = run.arm_config('fast')
    assert fast['phased_local_path'] + fast['phased_vcf_file_name'] == '../input/vcf/chr{chr}/*'


def test_subset_metadata_copies_lines(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    one_vcf(tmp_path)
    text = ('Sample name\tSex\tPopulation code\tSuperpopulation code\tScore\n'
            'S01\tmale\tFIN\tEUR\t1.0\nS02\tfemale\tGBR\tEUR\t \nS03\tfemale\tYRI\tAFR\t0012\n\n')
    (tmp_path / 'samples.tsv').write_text(text)
    run = verify.Run(config(), args(), str(tmp_path))
    os.makedirs(run.input_dir)
    assert verify.subset_metadata(run, ['S03', 'S01', 'S99']) == 2
    assert (Path(run.input_dir) / 'samples.tsv').read_text() == \
        'Sample name\tSex\tPopulation code\tSuperpopulation code\tScore\nS01\tmale\tFIN\tEUR\t1.0\n' \
        'S03\tfemale\tYRI\tAFR\t0012\n'
    assert run.arm_config('fast')['sample_metadata_file'] == '../input/samples.tsv'


def test_subset_transcripts_first_n_per_chromosome(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'tr.csv').write_text('chromosome,transcriptID\n' + ''.join(
        '%s,ENST%d\n' % (c, i) for i, c in enumerate([21, 22, 21, 21, 'X', 22, 21, 1])))
    run = verify.Run(config(custom_transcript_list='tr.csv', chromosomes=['21', 22, 'X']), args(),
                     str(tmp_path))
    os.makedirs(run.input_dir)
    assert verify.subset_transcripts(run, 2, None) == {'21': 2, '22': 2, 'X': 1}
    assert (Path(run.input_dir) / 'transcripts_first2.csv').read_text() == \
        'chromosome,transcriptID\n21,ENST0\n22,ENST1\n21,ENST2\nX,ENST4\n22,ENST5\n'
    assert run.arm_config('fast')['custom_transcript_list'] == '../input/transcripts_first2.csv'


# ------------------------------------ preflight ------------------------------------

def check(tmp_path, chunks, names=None):
    by_chrom = {'22': [str(write(tmp_path / n, t)) for n, t in chunks]}
    return preflight.check(by_chrom, names)


def test_preflight_clean_input(tmp_path):
    w, n = check(tmp_path, [('a.vcf', vcf_text(['A', 'B'], [(100, 'rs1', ['0|1', '1']), (200, 'rs2', ['.|0', '1|1'])]))],
                 {'A', 'B'})
    assert w == [] and n == []


@pytest.mark.parametrize('record,message', [
    ('22\t100\trs1\tA\tG\t.\tPASS\t.\tGT:DS\t0|1:1\t0|0:0\n', 'FORMAT fields other than GT'),
    ('22\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\t0|0\n', 'unphased'),
    ('22\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t10|1\t0|0\n', 'more than one character'),
    ('22\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\r\n', 'characters other than'),
    ('22\t100\trs1\tA\tG\t.\tPASS\t.\tDS\t1\t0\n', 'no GT in FORMAT'),
])
def test_preflight_genotype_warnings(tmp_path, record, message):
    w, _ = check(tmp_path, [('a.vcf', vcf_text(['A', 'B'], [record]))])
    assert len(w) == 1 and message in w[0]


def test_preflight_samples(tmp_path):
    w, _ = check(tmp_path, [('a.vcf', vcf_text(['123', 'B', 'B'], [(100, 'rs1', ['0|1'] * 3)])),
                            ('b.vcf', vcf_text(['A', 'B', 'C'], [(200, 'rs2', ['0|1'] * 3)]))], {'B', 'C'})
    text = '\n'.join(w)
    assert 'more than once: B' in text
    assert 'sample names are numbers (123)' in text
    assert 'lists other samples than' in text
    assert '1 of 3 samples are not in sample_metadata_file (123)' in text
    assert '1 of 3 samples are not in sample_metadata_file (A)' in text


def test_preflight_dot_ids_and_order(tmp_path):
    w, n = check(tmp_path, [('a.vcf', vcf_text(['A'], [(500, '.', ['0|1']), (400, 'rs', ['0|1'])])),
                            ('b.vcf', vcf_text(['A'], [(450, '.', ['0|1'])]))])
    text = '\n'.join(w)
    assert 'a.vcf is not sorted by position' in text
    assert 'a.vcf sorts before' in text and 'b.vcf' in text
    assert len(n) == 1 and 'records with ID "."' in n[0] and '(1 of 2)' in n[0] and '(1 of 1)' in n[0]


def test_preflight_overlap(tmp_path, monkeypatch):
    w, _ = check(tmp_path, [('a.vcf', vcf_text(['A'], [(100, 'rs1', ['0|1']), (300, 'rs2', ['0|1'])])),
                            ('b.vcf', vcf_text(['A'], [(200, 'rs3', ['0|1'])]))])
    assert len(w) == 1 and 'overlaps' in w[0] and 'ends at 300' in w[0]
    # a chunk scanned only in part: the index decides, or a note says the check needs one
    monkeypatch.setattr(preflight, 'MAX_RECORDS', 1)
    w, n = check(tmp_path, [('a.vcf', vcf_text(['A'], [(100, 'rs1', ['0|1']), (300, 'rs2', ['0|1'])])),
                            ('b.vcf', vcf_text(['A'], [(200, 'rs3', ['0|1'])]))])
    assert w == [] and len(n) == 1 and 'could not check whether' in n[0]
    tabix = shutil.which('tabix', path=str(BIN) + os.pathsep + os.environ.get('PATH', ''))
    bgzip = shutil.which('bgzip', path=str(BIN) + os.pathsep + os.environ.get('PATH', ''))
    if tabix and bgzip:
        for name, recs in (('c.vcf.gz', [(100, 'rs1', ['0|1']), (300, 'rs2', ['0|1'])]),
                           ('d.vcf.gz', [(200, 'rs3', ['0|1'])])):
            (tmp_path / name[:-3]).write_text(vcf_text(['A'], recs))
            subprocess.run([bgzip, '-f', str(tmp_path / name[:-3])], check=True)
            subprocess.run([tabix, '-p', 'vcf', str(tmp_path / name)], check=True)
        w, n = preflight.check({'22': [str(tmp_path / 'c.vcf.gz'), str(tmp_path / 'd.vcf.gz')]}, None, tabix)
        assert len(w) == 1 and 'has records at or after position 200' in w[0]


# ------------------------------------ dry run ------------------------------------

def test_dry_run_against_the_tiny_release(tmp_path, monkeypatch):
    pytest.importorskip('gffutils')
    if not (BIN / 'snakemake').exists():
        pytest.skip('snakemake is not installed next to ' + sys.executable)
    import tiny_ensembl as tiny
    downloads = tiny.write_downloads(tmp_path / 'downloads')
    proc = subprocess.run([sys.executable, str(ROOT / 'container' / 'prebuild_refs.py'), '--release', tiny.RELEASE,
                           '--downloads', str(downloads), '--out', str(tmp_path / 'refs'),
                           '--chromosomes', ','.join(tiny.CHROMOSOMES), '-j', '2'], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    work = tmp_path / 'work'
    shutil.copy(ROOT / 'crap.fasta', tmp_path / 'crap.fasta')
    write(work / 'vcfs' / 'chr21.vcf', vcf_text(SAMPLES, [(1100, '.', ['0|1'] * 40)]))
    write(work / 'vcfs' / 'chr22.vcf', vcf_text(SAMPLES, [(5100, '.', ['0|1'] * 40)]))
    (work / 'samples.tsv').write_text('Sample name\tSex\tPopulation code\tSuperpopulation code\n' + ''.join(
        '%s\t%s\tFIN\tEUR\n' % (s, 'male' if i % 2 else 'female') for i, s in enumerate(SAMPLES)))
    cfg = config(refs_dir=str(tmp_path / 'refs'), contaminants_fasta='../crap.fasta', included_transcript_biotypes='',
                 fasta_simplify_headers=True, working_dir_name_haplo='haplotypes_tmp', phased_min_af=0,
                 phased_af_field='AF', haplo_min_freq=0, haplo_min_count=0, x_par1_to=2781479,
                 x_par2_from=155701383, haplo_require_start=1, haplo_ignore_UTR=1, haplo_skip_start_lost=1,
                 max_cores=1, working_dir_name_var='variants_tmp', var_require_start=1, add_existing_haplo=False)
    with open(work / 'config.yaml', 'w') as fh:
        yaml.safe_dump(cfg, fh)
    monkeypatch.chdir(work)
    monkeypatch.delenv('PROHAP_REFS', raising=False)
    assert verify.main(['--cores', '2', '--chromosomes', '22', '--samples', '9', '--transcripts', '1', '--dry-run']) == 0
    run_dir = next((work / 'verify_runs').iterdir())
    fast = yaml.safe_load((run_dir / 'fast' / 'config.yaml').read_text())
    assert fast['chromosomes'] == ['22']
    assert (run_dir / 'input' / 'transcripts_first1.csv').read_text() == 'chromosome,transcriptID\n22,ENST00000000221\n'
    assert 'filter_phased_vcf' in (run_dir / 'fast.dryrun.log').read_text()


def test_vcfcut_changes_nothing_else(tmp_path):
    """Headers without definitions, undefined INFO tags, odd floats and CR line ends survive; bcftools view -s
    rewrites the first three and refuses undefined tags."""
    bare = '##fileformat=VCFv4.2\n##comment=a\ttab\n'
    text = vcf_text(['A', 'B', 'C'], ['22\t100\t.\tA\tG,T\t.\tPASS\tAF=1e-05,0.30000000000000004;XX=y\tGT\t0|1\t2|0\t1\r\n',
                                      '22\t200\trs2\tC\tA\t50\t.\t.\tGT\t.|.\t0|0\t1|1\n'], header=bare)
    write(tmp_path / 'in.vcf.gz', text)
    assert vcfcut.cut(str(tmp_path / 'in.vcf.gz'), str(tmp_path / 'out.vcf'), {b'A', b'C'}) == 2
    expected = text.replace('FORMAT\tA\tB\tC', 'FORMAT\tA\tC') \
        .replace('\t0|1\t2|0\t1\r\n', '\t0|1\t1\r\n').replace('\t.|.\t0|0\t1|1\n', '\t.|.\t1|1\n')
    assert (tmp_path / 'out.vcf').read_bytes() == expected.encode()


def test_vcfcut_errors(tmp_path):
    write(tmp_path / 'short.vcf', vcf_text(['A', 'B'], ['22\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\n']))
    with pytest.raises(ValueError, match='a line has 10 columns where the header has 11'):
        vcfcut.cut(str(tmp_path / 'short.vcf'), str(tmp_path / 'out.vcf'), {b'A'})
    with pytest.raises(ValueError, match='none of the selected samples'):
        vcfcut.cut(str(tmp_path / 'short.vcf'), str(tmp_path / 'out.vcf'), {b'Z'})
    write(tmp_path / 'nohead.vcf', '##fileformat=VCFv4.2\n22\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\n')
    with pytest.raises(ValueError, match='no #CHROM header line'):
        vcfcut.cut(str(tmp_path / 'nohead.vcf'), str(tmp_path / 'out.vcf'), {b'A'})
