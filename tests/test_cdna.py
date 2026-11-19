"""src/modules/cdna.py against the reference's read_fasta(filename, True):

    pytest tests/test_cdna.py
"""
import gzip
import os
import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))
from modules.cdna import read_sequences  # noqa: E402
from modules.legacy_import import legacy  # noqa: E402

REFS = Path(os.environ.get('PROHAP_REFS', '/nonexistent')).expanduser() / 'ensembl.v113'


def check(path):
    expected = {k: v['sequence'] for k, v in legacy.common.read_fasta(str(path), True).items()}
    assert read_sequences(str(path), set(expected) | {'absent'}) == expected
    some = set(list(expected)[::2])
    assert read_sequences(str(path), some) == {k: expected[k] for k in some}
    return expected


@pytest.mark.parametrize('ending', ['\n', '', '\n\n'])
@pytest.mark.parametrize('gz', [False, True])
def test_crafted(tmp_path, ending, gz):
    rng = random.Random(len(ending) + gz)
    headers = ['>ENST0001.2 cdna chromosome:GRCh38:22:1:100:1 gene:ENSG1', '>ENST0002 cdna', '>generic_x|ENST0003.1|desc',
               '>sp|ENST0004.5', '>sp|ENST0005.1 with space', '>a|b|c|d', '>ENST0001.9 duplicate, replaces the first',
               '>', '>ENST0006.1\tsecond']
    lines = []
    for h in headers:
        lines.append(h)
        for _ in range(rng.randint(0, 3)):
            lines.append(''.join(rng.choice('ACGTN') for _ in range(rng.randint(0, 70))))
    text = '\n'.join(lines) + ending
    path = tmp_path / ('c.fa.gz' if gz else 'c.fa')
    if gz:
        with gzip.open(path, 'wt') as f:
            f.write(text)
    else:
        path.write_text(text)
    expected = check(path)
    assert 'ENST0001' in expected and len(expected) >= 7


def test_odd_files(tmp_path):
    for text in ['', 'ACGT\n>T1\nAC\n', '\n>T1\nAC\n', '>T1', '>T1\n', '>T1\nA', '>T1\r\nACG\r\nT\r\n>T2\rAC\r']:
        path = tmp_path / 'odd.fa'
        path.write_bytes(text.encode())
        check(path)


@pytest.mark.skipif(not REFS.is_dir(), reason='needs the v113 references in $PROHAP_REFS')
def test_ensembl_cdna():
    expected = check(REFS / 'total_cdnas_113.fa')
    assert len(expected) > 200000
