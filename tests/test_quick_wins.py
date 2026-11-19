"""The quick wins in src/ (IMPLEMENTATION_NOTES section 1) against legacy/src, function by function, on small
crafted inputs. tests/legacy_vs_fast.sh is the real check; these run in the unit CI, in a few seconds:

    pytest tests/test_quick_wins.py
"""
import importlib
import io
import os
import random
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import pytest

pd = pytest.importorskip('pandas')

ROOT = Path(__file__).resolve().parent.parent
TREES = {'legacy': ROOT / 'legacy' / 'src', 'fast': ROOT / 'src'}
VCF_COLUMNS = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
PAR1_TO, PAR2_FROM = 2781479, 155701383


@contextmanager
def modules_of(arm):
    """Import `modules.<name>` from one tree. `modules` is a namespace package in both trees, so only one tree may be
    on sys.path at a time and its imports are dropped afterwards."""
    def drop():
        for name in [n for n in sys.modules if n == 'modules' or n.startswith('modules.')]:
            del sys.modules[name]
    drop()
    sys.path.insert(0, str(TREES[arm]))
    try:
        yield lambda name: importlib.import_module('modules.' + name)
    finally:
        sys.path.remove(str(TREES[arm]))
        drop()


class SerialPool:
    """Stands in for multiprocessing.Pool: the reference's closures need fork, which macOS does not use by default."""
    def __init__(self, processes=None):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def map(self, func, items):
        return [func(item) for item in items]


def transcript_table(path, rows, samples):
    """A per-transcript table as vcf_reader.parse_vcf writes it: the VCF lines read by pandas and written back."""
    lines = ['\t'.join(VCF_COLUMNS + samples)]
    for pos, vcf_id, ref, alt, info, genotypes in rows:
        lines.append('\t'.join(['22', str(pos), vcf_id, ref, alt, '.', 'PASS', info, 'GT'] + genotypes))
    df = pd.read_csv(io.StringIO('\n'.join(lines) + '\n'), sep='\t')
    df.to_csv(path, sep='\t', index=False, header=True)


# ---------------------------------------------------------------------------------------------------------------
# P1: get_haplotypes.py reads POS/REF/ALT/ID/INFO as lists once instead of vcf_df.iloc[idx] per haplotype and variant


def wide_table(ids):
    """A table 62,009 columns wide: 2**20 // 62009 = 16, so pandas types it in chunks of 8 rows, as at 100k samples."""
    samples = '\t'.join('0|0' for _ in range(62000))
    lines = ['\t'.join(VCF_COLUMNS + ['S%d' % i for i in range(62000)])]
    for i, vcf_id in enumerate(ids):
        lines.append('\t'.join(['22', str(100 + i), vcf_id, 'A', 'G', '.', 'PASS', 'AC=1;AN=2;AF=0.5', 'GT', samples]))
    return pd.read_csv(io.StringIO('\n'.join(lines) + '\n'), sep='\t')


@pytest.mark.filterwarnings('ignore::pandas.errors.DtypeWarning')
@pytest.mark.parametrize('ids, dtype', [
    # chunks of ints, floats, booleans and strings: one object column holding all four Python types
    (['7', '123', '0012', '5', '6', '8', '9', '10',
      '2e5', '1e3', '7', 'NA', 'nan', '1.50', '3e4', '12',
      'True', 'False', 'TRUE', 'true', 'False', 'True', 'FALSE', 'false',
      'rs1', '2e5', 'NA', 'True', '0012', '7', 'abc', '1a2b'], 'object'),
    # chunks of ints and floats: a float64 column, numpy scalars from iloc
    (['7', '123', '0012', '5', '6', '8', '9', '10',
      '2e5', '1e3', '1e16', 'NA', 'nan', '1.50', '3e-5', '12'], 'float64'),
])
def test_p1_list_values_print_like_iloc_values(ids, dtype):
    df = wide_table(ids)
    assert str(df['ID'].dtype) == dtype
    columns = {c: df[c].tolist() for c in ('POS', 'REF', 'ALT', 'ID', 'INFO')}
    for i in range(len(df)):
        row = df.iloc[i]
        for c, values in columns.items():
            assert str(values[i]) == str(row[c]), (i, c)


def haplotype_inputs(tmp_path):
    """Transcripts, per-transcript tables and sample metadata that reach every branch of get_gene_haplotypes: IDs that
    pandas types as int64, float64, bool and object, INFO without AF, a deletion overlapping an SNV on one copy
    (conflict, removed samples), a transcript without variants, PAR1 and non-PAR transcripts on chromosome X, a VCF
    sample missing from the metadata and a metadata row missing from the VCF."""
    rng = random.Random(7)
    samples = ['S%02d' % i for i in range(1, 13)]
    metadata = pd.DataFrame({
        'Sample name': samples[:-1] + ['M01'],
        'Sex': ['male' if i % 2 else 'female' for i in range(1, 13)],
        'Population code': [rng.choice(['GBR', 'YRI', 'CHB']) for _ in range(12)],
        'Superpopulation code': '',
    })
    metadata['Superpopulation code'] = metadata['Population code'].map({'GBR': 'EUR', 'YRI': 'AFR', 'CHB': 'EAS'})
    indiv_ids = [s for s in samples if s in set(metadata['Sample name'])]

    def genotypes(haploid_males=False):
        tokens = []
        for i, _ in enumerate(samples, start=1):
            a, b = rng.random() < 0.3, rng.random() < 0.3
            tokens.append('%d|%d' % (a, 0 if (haploid_males and i % 2) else b))
        return tokens

    tables = {
        # object column: rs IDs with tokens pandas would retype on their own
        'T1': [(1100, '2e5', 'A', 'G', 'AC=3;AN=24;AF=0.125', genotypes()),
               (1200, 'NA', 'C', 'T', 'AC=1;AN=24;AF=0.04', genotypes()),
               (1300, 'rs7', 'GAT', 'G', 'AC=2;AN=24;AF=0.08', genotypes()),
               (1301, '0012', 'A', 'C', 'AC=2;AN=24', genotypes()),
               (1400, 'rs9', 'T', 'TA', 'AC=2;AN=24;AF=0.08;DP=7', genotypes())],
        'T2': [],
        'T3': [(2100, 'True', 'A', 'G', 'AF=0.2', genotypes()),      # bool column
               (2200, 'False', 'C', 'G', 'AF=0.3', genotypes())],
        'T4': [(3100, '7', 'A', 'G', 'AF=0.2', genotypes()),         # int64 column
               (3200, '0012', 'C', 'G', 'AF=0.3', genotypes())],
        'T5': [(3000100, '1e3', 'A', 'G', 'AF=0.2', genotypes(True)),     # float64 column, beyond PAR1
               (3000200, 'nan', 'C', 'G', 'AF=0.3', genotypes(True)),
               (3000300, '7', 'G', 'T', 'AF=0.1', genotypes(True))],
    }
    spans = {'T1': (1000, 5000), 'T2': (6000, 7000), 'T3': (2000, 2500), 'T4': (3000, 3500),
             'T5': (3000000, 3001000)}
    for tid, rows in tables.items():
        transcript_table(tmp_path / (tid + '.tsv'), rows, samples)
    transcripts = [SimpleNamespace(id=tid, start=s, end=e) for tid, (s, e) in spans.items()]
    return transcripts, indiv_ids, metadata


def run_get_gene_haplotypes(arm, tmp_path, monkeypatch, is_x, edit=None):
    transcripts, indiv_ids, metadata = haplotype_inputs(tmp_path)
    if edit:
        indiv_ids, metadata = edit(indiv_ids, metadata)
    log = tmp_path / (arm + '.log')
    with modules_of(arm) as module:
        get_haplotypes = module('get_haplotypes')
        monkeypatch.setattr(get_haplotypes, 'Pool', SerialPool)
        df = get_haplotypes.get_gene_haplotypes(transcripts, indiv_ids, str(tmp_path), str(log), 1, is_x, PAR1_TO,
                                                PAR2_FROM, metadata.copy())
    return df.to_csv(sep='\t', index=False), log.read_text()


@pytest.mark.parametrize('is_x', [False, True], ids=['autosome', 'chrX'])
def test_get_gene_haplotypes_matches_legacy(tmp_path, monkeypatch, is_x):
    legacy_table, legacy_log = run_get_gene_haplotypes('legacy', tmp_path, monkeypatch, is_x)
    fast_table, fast_log = run_get_gene_haplotypes('fast', tmp_path, monkeypatch, is_x)
    assert fast_table == legacy_table
    assert fast_log == legacy_log
    # the inputs reach what they are meant to: IDs as pandas typed them, conflicting variants
    ids = {i for line in legacy_table.splitlines()[1:] for i in line.split('\t')[3].split(';')}
    assert {'2e5', '0012', 'nan', 'True', '12', '1000.0', '7.0'} <= ids
    assert any(line.split('\t')[1] != '0' for line in legacy_log.splitlines()[1:])


# ---------------------------------------------------------------------------------------------------------------
# P2: get_haplotypes.py looks up sex and populations in dicts instead of sample_info.loc per sample


@pytest.mark.parametrize('is_x', [False, True], ids=['autosome', 'chrX'])
def test_vcf_sample_listed_twice(tmp_path, monkeypatch, is_x):
    """prohap.py takes the sample names from the VCF header as they are, so a name listed twice is counted twice (and
    pandas reads the first of its columns both times)."""
    def edit(indiv_ids, metadata):
        return indiv_ids + [indiv_ids[0]], metadata
    legacy = run_get_gene_haplotypes('legacy', tmp_path, monkeypatch, is_x, edit)
    assert run_get_gene_haplotypes('fast', tmp_path, monkeypatch, is_x, edit) == legacy


def test_metadata_name_listed_twice(tmp_path, monkeypatch):
    """The reference stops when a sample of the VCF has two metadata rows; so must src/."""
    def edit(indiv_ids, metadata):
        return indiv_ids, pd.concat([metadata, metadata.iloc[[2]]], ignore_index=True)
    for arm in ('legacy', 'fast'):
        with pytest.raises(ValueError, match='truth value of a Series is ambiguous'):
            run_get_gene_haplotypes(arm, tmp_path, monkeypatch, True, edit)


# ---------------------------------------------------------------------------------------------------------------
# P4: vcf_reader.py collects each transcript's VCF lines in a list and joins them once


class FakeDB:
    """The one gffutils call parse_vcf makes: the exons of a transcript, ordered by start."""
    def __init__(self, exons):
        self.exons = exons

    def children(self, transcript, featuretype, order_by):
        assert (featuretype, order_by) == ('exon', 'start')
        return [SimpleNamespace(start=s, end=e) for s, e in self.exons[transcript.id]]


def test_parse_vcf_matches_legacy(tmp_path):
    rng = random.Random(3)
    exons = {'T1': [(1000, 1100), (1500, 1600)], 'T2': [(1050, 1080), (1550, 1700)], 'T3': [(1520, 1530)],
             'T4': [(5000, 5100)], 'T5': [(9000, 9100)]}
    transcripts = [SimpleNamespace(id=t, start=e[0][0], end=e[-1][1]) for t, e in exons.items()]
    samples = ['S%d' % i for i in range(1, 9)]
    lines = ['##fileformat=VCFv4.2', '##contig=<ID=22>', '#' + '\t'.join(VCF_COLUMNS + samples)]
    positions = [500, 1000, 1001, 1060, 1098, 1100, 1300, 1495, 1525, 1525, 1560, 1650, 3000, 5050, 5100, 7000]
    for n, pos in enumerate(positions):
        ref, alt = rng.choice([('A', 'G'), ('C', 'T'), ('ACGTA', 'A'), ('G', 'GTT'), ('T', 'C,G'), ('N', 'A')])
        vcf_id = rng.choice(['.', '.', 'rs%d' % n, '2e5', 'NA'])
        af = rng.choice(['0.01', '0.2', '0.5'])
        gts = ['%d|%d' % (rng.random() < 0.3, rng.random() < 0.3) for _ in samples]
        lines.append('\t'.join(['22', str(pos), vcf_id, ref, alt, '.', 'PASS', 'AC=1;AN=16;AF=' + af, 'GT'] + gts))
    text = '\n'.join(lines) + '\n'

    written = {}
    for arm in ('legacy', 'fast'):
        out = tmp_path / arm
        out.mkdir()
        with modules_of(arm) as module:
            colnames = module('vcf_reader').parse_vcf(transcripts, io.StringIO(text), FakeDB(exons), 0.05, str(out))
        written[arm] = colnames, {p.name: p.read_bytes() for p in sorted(out.iterdir())}
    assert written['fast'] == written['legacy']
    tables = written['legacy'][1]
    assert sorted(tables) == ['T%d.tsv' % i for i in range(1, 6)]
    assert sum(len(t.splitlines()) - 1 for t in tables.values()) > len(exons)    # lines shared by transcripts


# ---------------------------------------------------------------------------------------------------------------
# D3: merge_duplicate_seq.py merges by sequence hash in a dict and sorts once


def test_merge_duplicate_seq_matches_legacy(tmp_path):
    """Both scripts on 60 records of 15 sequences, with every kind of header the script parses. The record order is
    by hash(), so the scripts run with PYTHONHASHSEED=0 as in the workflow."""
    rng = random.Random(5)
    sequences = [''.join(rng.choice('ACDEFGHIKLMN') for _ in range(rng.randint(5, 30))) for _ in range(15)]
    headers = [
        'generic_ensref|ENSP{i:05d}|matching_proteins:ENST{i:05d} position_within_protein:{r}',
        'generic_enshap|enshap_22_{i:x}|matching_proteins:haplo_chr22_{i:x};haplo_chr22_{j:x} start:{r} '
        'reading_frame:{r};{s}',
        'generic_var|var_{i}|matching_proteins:var_{i} start:0 reading_frame:{r} position_within_protein:{s}',
        'generic_cont|CONT{i}|contaminant {i}',
        'sp|P{i:05d}|Swiss-Prot entry {i}',
    ]
    fasta = tmp_path / 'input.fa'
    with open(fasta, 'w') as out:
        for i in range(60):
            header = headers[i % 5].format(i=i, j=i + 100, r=rng.randint(0, 2), s=rng.randint(0, 2))
            out.write('>' + header + '\n' + rng.choice(sequences) + '\n')

    env = dict(os.environ, PYTHONHASHSEED='0', PYTHONNOUSERSITE='1')
    outputs = {}
    for arm in ('legacy', 'fast'):
        result = tmp_path / (arm + '.fa')
        subprocess.run([sys.executable, str(TREES[arm] / 'merge_duplicate_seq.py'), '-i', str(fasta), '-o',
                        str(result)], env=env, check=True)
        outputs[arm] = result.read_text()
    assert outputs['fast'] == outputs['legacy']
    assert outputs['legacy'].count('>') == len(set(sequences))
