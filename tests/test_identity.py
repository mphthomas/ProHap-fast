"""Identity tests: verify (legacy arm vs fast arm) must find every output IDENTICAL on the fixtures, and DIFFERENT
when the fast arm's code changes one character of an output.

    PROHAP_REFS=~/prohap-refs pytest tests/test_identity.py

About 3 minutes in all on 4 cores. tests/legacy_vs_fast.sh is the check of src/ against legacy/src, on these
fixtures and simulated ones; this module tests verify, the user-facing check.
Needs the ProHap environment and the pre-built references of release 113 ($PROHAP_REFS/ensembl.v113, else
/opt/prohap/refs/ensembl.v113); skipped without them. VERIFY_CORES sets the cores for each run (default 4).

The fixtures are described in tests/fixtures/1kg/README.md.
"""
import os
import re
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parent.parent
FIXTURES = ROOT / 'tests' / 'fixtures' / '1kg'
TARBALL = ROOT / 'sample_data' / 'sample_1kGP_common_global.tar.gz'
VERIFY = ROOT / 'verify' / 'verify.py'
CORES = os.environ.get('VERIFY_CORES', '4')
REFS = Path(os.environ.get('PROHAP_REFS') or '/opt/prohap/refs').expanduser()

if not (Path(sys.executable).parent / 'snakemake').exists():
    pytest.skip('snakemake is not installed next to ' + sys.executable, allow_module_level=True)
if not (REFS / 'ensembl.v113' / 'release.txt').exists():
    pytest.skip('no pre-built references for release 113 in ' + str(REFS) + ' (set PROHAP_REFS)',
                allow_module_level=True)


@pytest.fixture(scope='module')
def sample_data(tmp_path_factory):
    """The 1kG VCFs of chromosomes 1 (two chunks), 21, 22 and X, extracted once."""
    dest = tmp_path_factory.mktemp('sample_data')
    wanted = re.compile(r'sample_chr(1|21|22|X)_1kGP_GRCh38_phased(_2)?\.vcf')
    with tarfile.open(TARBALL) as tar:
        members = [m for m in tar.getmembers() if m.isfile() and wanted.fullmatch(m.name)]
        assert len(members) == 5
        tar.extractall(dest, members=members, filter='data')
    return dest


def workdir(tmp_path, sample_data, fixture, **overrides):
    """A work directory as a user would set it up: config.yaml, the VCFs, the metadata and the contaminants."""
    work = tmp_path / 'work'
    work.mkdir()
    (work / 'sample_data').symlink_to(sample_data)
    (work / 'igsr_samples.tsv').symlink_to(ROOT / 'igsr_samples.tsv')
    (work / 'crap.fasta').symlink_to(ROOT / 'crap.fasta')
    (work / 'transcripts.csv').symlink_to(FIXTURES / 'transcripts.csv')
    text = (FIXTURES / (fixture + '.yaml')).read_text()
    if overrides:
        cfg = yaml.safe_load(text)
        cfg.update(overrides)
        text = yaml.safe_dump(cfg, sort_keys=False)
    (work / 'config.yaml').write_text(text)
    return work


def verify(work, *args):
    env = dict(os.environ, PROHAP_REFS=str(REFS))
    proc = subprocess.run([sys.executable, str(VERIFY), '--cores', CORES] + list(args), cwd=work, env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(proc.stdout)
    return proc


def table(output):
    """{file: result} from verify's table."""
    rows = {}
    for line in output.splitlines():
        fields = line.split()
        if len(fields) >= 4 and fields[3] in ('IDENTICAL', 'DIFFERENT', 'MISSING'):
            rows[fields[0]] = fields[3]
    return rows


@pytest.mark.parametrize('fixture', ['threshold0', 'template'])
def test_fixture_is_identical(tmp_path, sample_data, fixture):
    work = workdir(tmp_path, sample_data, fixture)
    proc = verify(work)
    assert proc.returncode == 0
    rows = table(proc.stdout)
    # 4 per chromosome (filtered VCF, table, FASTA, log) and 11 others, plus the cDNA FASTA with threshold0
    assert len(rows) == 4 * 4 + 11 + (fixture == 'threshold0')
    assert set(rows.values()) == {'IDENTICAL'}
    assert 'IDENTICAL: all %d outputs' % len(rows) in proc.stdout
    assert 'the legacy arm runs at full scale' in proc.stdout
    # only the report is left
    assert [p.name for p in (work / 'verify_runs').glob('*/*')] == ['report.tsv']


def test_subsets_are_identical(tmp_path, sample_data):
    """Random samples and the first transcripts of chromosome 1, which comes in two chunks, and 22."""
    work = workdir(tmp_path, sample_data, 'threshold0', chromosomes=['1', '22'])
    proc = verify(work, '--samples', '300', '--transcripts', '30', '--keep')
    assert proc.returncode == 0
    assert set(table(proc.stdout).values()) == {'IDENTICAL'}
    run = next((work / 'verify_runs').iterdir())
    assert sorted(os.listdir(run / 'input' / 'vcf' / 'chr1')) == ['sample_chr1_1kGP_GRCh38_phased.vcf',
                                                                   'sample_chr1_1kGP_GRCh38_phased_2.vcf']
    assert len((run / 'input' / 'samples.txt').read_text().split()) == 300
    assert len((run / 'input' / 'igsr_samples.tsv').read_text().splitlines()) == 301
    assert (run / 'legacy' / 'results' / '1kg_haplotypes.tsv').stat().st_size > 1000


def test_changed_fast_code_is_different(tmp_path, sample_data):
    """A copy of src/ that writes one character differently in the removed-samples log: verify must report that
    file DIFFERENT, everything else IDENTICAL, and exit with status 1."""
    fast_src = tmp_path / 'src'
    shutil.copytree(ROOT / 'src', fast_src, ignore=shutil.ignore_patterns('__pycache__'))
    module = fast_src / 'modules' / 'get_haplotypes.py'
    code = module.read_text()
    assert code.count("'TranscriptID\\t#Removed\\tRemovedSamples\\n'") == 1
    module.write_text(code.replace("'TranscriptID\\t#Removed", "'TranscriptId\\t#Removed"))
    work = workdir(tmp_path, sample_data, 'threshold0', chromosomes=['22'])
    proc = verify(work, '--transcripts', '20', '--fast-src', str(fast_src))
    assert proc.returncode == 1
    rows = table(proc.stdout)
    assert rows.pop('log/prohap_chr22.log') == 'DIFFERENT'
    assert rows and set(rows.values()) == {'IDENTICAL'}
    assert 'log/prohap_chr22.log, line 1:\n  legacy: TranscriptID\t#Removed' in proc.stdout
    assert 'DIFFERENT: 1 of' in proc.stdout
    # everything is kept for inspection
    run = next((work / 'verify_runs').iterdir())
    assert (run / 'fast' / 'log' / 'prohap_chr22.log').read_text().startswith('TranscriptId')
