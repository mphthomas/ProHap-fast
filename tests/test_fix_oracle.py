"""The naive oracle (tests/fix_oracle.py) against prohap.py -fix_bugs 1, copy by copy, and the diff report's
comparison (tests/fix_diff_report.py), on the crafted chromosome of test_fix_bugs.py and, with the release 113
references, on chromosome 22 of the 1kG fixture. The whole set (the 1kG fixture and the simulated mix, TTN and
benchmark-block data, about 2 minutes):

    python3 tests/fix_oracle.py all

    pytest tests/test_fix_oracle.py
"""
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
pytest.importorskip('gffutils')
import fix_diff_report  # noqa: E402
import fix_oracle  # noqa: E402
import test_fix_bugs  # noqa: E402

REFS = Path(os.environ.get('PROHAP_REFS') or '/opt/prohap/refs').expanduser() / 'ensembl.v113'


def test_crafted_chromosome(tmp_path):
    """Every case of test_fix_bugs.py: the oracle's protein (or removal) is fast's, except the copy whose REF is not
    in the cDNA (left out: the oracle does not model it)."""
    test_fix_bugs.write_inputs(tmp_path)
    out = tmp_path / 'fast'
    test_fix_bugs.run(tmp_path, fix_oracle.ROOT / 'src' / 'prohap.py', out, 1)
    counts, shown = fix_oracle.compare(out, tmp_path / 'in.vcf.gz', '22', tmp_path / 'samples.tsv',
                                       tmp_path / 'transcripts.csv', db_path=tmp_path / 'in.db',
                                       cdna_path=tmp_path / 'cdna.fa')
    assert counts['DIFFERENT'] == 0 and counts['shared_id'] == 0, shown
    assert counts['left out: REF not in the cDNA'] == 1
    # every case's copy and the reference copies of both transcripts
    assert counts['same'] == counts['copies compared'] == 2 * 2 * len(test_fix_bugs.NAMES) - 1


def test_diff_report(tmp_path):
    """The diff report's counts between the reference's output and the fixed output of the crafted chromosome."""
    test_fix_bugs.write_inputs(tmp_path)
    test_fix_bugs.run(tmp_path, fix_oracle.ROOT / 'src' / 'prohap.py', tmp_path / 'off', 0)
    test_fix_bugs.run(tmp_path, fix_oracle.ROOT / 'src' / 'prohap.py', tmp_path / 'on', 1)
    counts = fix_diff_report.diff(fix_diff_report.read_run(tmp_path / 'off'), fix_diff_report.read_run(tmp_path / 'on'))
    removed = [n for n in test_fix_bugs.NAMES if test_fix_bugs.CASES[n][2]['off'] == 'removed']
    stay = [n for n in removed if test_fix_bugs.expected(n, test_fix_bugs.fixes.ALL) == 'removed']
    assert counts['removed copies: no longer removed'] == len(removed) - len(stay)
    assert counts['removed copies: newly removed'] == 0
    assert counts['copies: removed -> protein'] == len(removed) - len(stay)
    assert counts['copies: protein -> other protein'] == 1          # P11c: the reference missed its lost start
    assert counts['copies: ref -> protein'] == sum(
        1 for n in test_fix_bugs.NAMES if test_fix_bugs.CASES[n][2]['off'] == 'ref' and
        test_fix_bugs.expected(n, test_fix_bugs.fixes.ALL) != 'ref')


@pytest.mark.skipif(not (REFS / 'release.txt').exists(), reason='needs the release 113 references (PROHAP_REFS)')
def test_1kg_chromosome_22(tmp_path):
    datasets = fix_oracle.datasets(tmp_path)
    chrom, vcf_glob, af, samples, transcripts = next(d for d in datasets['1kg']() if d[0] == '22')
    out = tmp_path / '1kg' / 'chr22' / 'fast'
    filtered = fix_oracle.run_prohap(out, chrom, vcf_glob, af, samples, transcripts)
    counts, shown = fix_oracle.compare(out, filtered, chrom, samples, transcripts)
    assert counts['DIFFERENT'] == 0, shown
    assert counts['same'] == counts['copies compared'] > 100000
