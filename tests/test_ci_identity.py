"""The CI identity plumbing (.github/workflows/identity.yml, tests/ci_identity.sh, legacy_vs_fast.sh -g/-G) and the
golden checksums, without an image: what can be checked here is that the pieces name each other and that the golden
file covers what CI checks.

    pytest tests/test_ci_identity.py
"""
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
GOLDEN = ROOT / 'tests' / 'fixtures' / 'expected.v113.sha256'
SCRIPT = ROOT / 'tests' / 'legacy_vs_fast.sh'


def tier(name):
    m = re.search(r'^%s="([^"]*)"' % name, SCRIPT.read_text(), re.M)
    return m.group(1).split()


def golden():
    entries = {}
    for line in GOLDEN.read_text().splitlines():
        if line.strip() and not line.startswith('#'):
            digest, name = line.split('  ', 1)
            assert re.fullmatch(r'[0-9a-f]{64}', digest), line
            entries[name] = digest
    return entries


def test_golden_covers_the_tiers_ci_runs():
    entries = golden()
    for fixture in tier('SMALL') + tier('MEDIUM'):
        names = [n for n in entries if n.startswith(fixture + '/')]
        assert names, 'no golden checksums for ' + fixture
        if fixture != 'dupseq':
            assert any(n.endswith('_simplified.fasta') for n in names), fixture
            assert not any('_phased_filtered.vcf.gz' in n for n in names), fixture
    assert not [n for n in entries if not any(n.startswith(f + '/') for f in tier('SMALL') + tier('MEDIUM')
                                              + tier('LARGE'))]


def test_workflow_runs_the_script():
    workflow = (ROOT / '.github' / 'workflows' / 'identity.yml').read_text()
    assert 'tests/ci_identity.sh' in workflow and (ROOT / 'tests' / 'ci_identity.sh').is_file()
    assert "'small medium'" in workflow and "'small'" in workflow
    assert 'apptainer' in workflow


@pytest.mark.skipif(sys.platform == 'linux', reason='the golden checksums are checked on Linux')
@pytest.mark.parametrize('flag', ['-g', '-G'])
def test_golden_needs_linux(flag):
    proc = subprocess.run(['bash', str(SCRIPT), flag, 'small'], capture_output=True, text=True)
    assert proc.returncode == 2 and 'need Linux' in proc.stderr
