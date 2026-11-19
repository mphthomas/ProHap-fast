"""src/modules/sweep.py (R1a) against the reference's parse_vcf: which records each transcript gets, in which order the
transcripts leave the queue, what get_haplotypes reads back from each transcript's table, and who carries what.

    pytest tests/test_sweep.py
"""
import random
from pathlib import Path
from types import SimpleNamespace

import pytest

pd = pytest.importorskip('pandas')

from fast_cases import COLUMNS, ROOT, fast, make_case, reference  # noqa: E402
from modules import sweep  # noqa: E402
from modules.fallback import Unsupported  # noqa: E402
from modules.pandas_emulation import pandas_chunk_rows  # noqa: E402


@pytest.mark.parametrize('seed, unsorted', [(1, False), (2, False), (3, True), (4, False), (5, True)])
def test_sweep_matches_parse_vcf(seed, unsorted, tmp_path, monkeypatch):
    transcripts, exons, samples, text = make_case(seed, unsorted=unsorted)
    if seed == 2:
        text = text.rstrip('\n')                       # no newline at the end of the last record
    included = [s for s in samples if s not in ('S03', 'S07')]          # two samples not in the metadata
    colnames, order, tables = reference(transcripts, exons, text, tmp_path, monkeypatch)
    columns, released = fast(transcripts, exons, samples, included, text, tmp_path)
    assert columns == colnames
    assert [transcripts[i].id for i, _ in released] == order
    chunk_rows = pandas_chunk_rows(len(columns))
    assigned = 0
    for i, records in released:
        df = tables[transcripts[i].id]
        assert len(records) == len(df)
        got = sweep.transcript_columns(records, chunk_rows)
        for column in ('POS', 'ID', 'REF', 'ALT', 'INFO'):
            assert got[column] == [str(v) for v in df[column].tolist()], (transcripts[i].id, column)
        # carriers: first copy where the token starts with '1|', second where it ends with '|1'
        expected = set()
        for rank, sample in enumerate(included):
            for row, token in enumerate(df[sample].tolist()):
                if token.startswith('1|'):
                    expected.add((row, 2 * rank))
                if token.endswith('|1'):
                    expected.add((row, 2 * rank + 1))
        assert {(row, int(c)) for row, r in enumerate(records) for c in r[6]} == expected
        assigned += len(records)
    assert assigned > 100


def test_empty_vcf(tmp_path):
    path = tmp_path / 'in.vcf'
    path.write_text('##fileformat=VCFv4.2\n#' + '\t'.join(COLUMNS[1:] + ['S1']) + '\n')
    reader = sweep.VCFReader(str(path))
    assert sweep.read_header(reader) == (3, None, b'')


@pytest.mark.parametrize('edit, message', [
    (lambda t: t.replace('\n', '\r\n'), 'carriage return'),
    (lambda t: t.replace('PASS', 'PÄSS', 1), 'non-ASCII'),
    (lambda t: t.replace('PASS', '"PASS"', 1), 'quote'),
    (lambda t: t.replace('\t0|1', '\t0/1', 1), r'a\|b'),
    (lambda t: t.replace('\t0|1', '\t10|1', 1), r'a\|b'),
    (lambda t: t.replace('\tGT\t', '\tGT:DS\t', 1), 'FORMAT'),
    (lambda t: t.replace('\t22\t', '\t022\t', 1).replace('\n22\t1', '\n22\t01', 1), 'POS'),
    (lambda t: t.replace('AC=1', 'AC=1 x', 1), 'whitespace'),
    (lambda t: t.replace('\tS01\t', '\tS00\t'), 'repeated'),
    (lambda t: t.replace('#CHROM', '#CHR'), 'nine VCF columns'),
    (lambda t: t.replace('AC=1;AN=24;AF=0.1', 'NA'), 'INFO'),
])
def test_unsupported_input(edit, message, tmp_path):
    rng = random.Random(9)
    samples = ['S%02d' % i for i in range(4)]
    lines = ['\t'.join(COLUMNS + samples)]
    for pos in range(100, 140, 2):
        lines.append('\t'.join(['22', str(pos), '.', 'A', 'G', '.', 'PASS', 'AC=1;AN=24;AF=0.1', 'GT']
                               + [rng.choice(['0|0', '0|1', '1|0']) for _ in samples]))
    lines[3] = lines[3].rsplit('\t', 1)[0] + '\t0|1'
    text = edit('\n'.join(lines) + '\n')
    transcripts = [SimpleNamespace(id='T1', start=90, end=200)]
    with pytest.raises(Unsupported, match=message):
        fast(transcripts, {'T1': [(90, 200)]}, samples, samples, text, tmp_path)


REFS = Path(__import__('os').environ.get('PROHAP_REFS', '/nonexistent')).expanduser() / 'ensembl.v113'


@pytest.mark.skipif(not REFS.is_dir(), reason='needs the v113 references in $PROHAP_REFS')
def test_sweep_on_1kg_chr22(tmp_path, monkeypatch):
    """The 1kG sample data of chromosome 22 (2,548 samples, `.` IDs) on the 1kG fixture's transcripts."""
    import tarfile
    gffutils = pytest.importorskip('gffutils')
    with tarfile.open(ROOT / 'sample_data' / 'sample_1kGP_common_global.tar.gz') as tar:
        text = tar.extractfile('sample_chr22_1kGP_GRCh38_phased.vcf').read().decode()
    ids = pd.read_csv(ROOT / 'tests' / 'fixtures' / '1kg' / 'transcripts.csv')
    ids = ids[ids['chromosome'].astype(str) == '22']['transcriptID'].tolist()
    db = gffutils.FeatureDB(str(REFS / 'Homo_sapiens.GRCh38.113.chr_patch_hapl_scaff_chr22.db'))
    features = sorted((db[t] for t in ids), key=lambda f: f.start)
    exons = {f.id: [(e.start, e.end) for e in db.children(f, featuretype='exon', order_by='start')]
             for f in features}
    transcripts = [SimpleNamespace(id=f.id, start=f.start, end=f.end) for f in features]
    samples = text.split('\n#CHROM', 1)[1].split('\n', 1)[0].split('\t')[9:]
    included = samples[::3]
    colnames, order, tables = reference(transcripts, exons, text, tmp_path, monkeypatch)
    columns, released = fast(transcripts, exons, samples, included, text, tmp_path)
    assert columns == colnames and [transcripts[i].id for i, _ in released] == order
    chunk_rows = pandas_chunk_rows(len(columns))
    assigned = 0
    for i, records in released:
        df = tables[transcripts[i].id]
        got = sweep.transcript_columns(records, chunk_rows)
        for column in ('POS', 'ID', 'REF', 'ALT', 'INFO'):
            assert got[column] == [str(v) for v in df[column].tolist()]
        for row, record in enumerate(records):
            tokens = df.iloc[row][included].tolist()
            expected = [2 * k for k, t in enumerate(tokens) if t.startswith('1|')] + \
                [2 * k + 1 for k, t in enumerate(tokens) if t.endswith('|1')]
            assert record[6].tolist() == sorted(expected)
        assigned += len(records)
    assert assigned > 30
