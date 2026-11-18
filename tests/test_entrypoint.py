"""The container entrypoint and container/check_release.py, run outside the image.

entrypoint.sh finds the image's files next to itself, so each test builds a stand-in image directory: the real
entrypoint.sh and check_release.py, a VERSION, refs/ensembl.v113/release.txt, and stubs for snakemake and
verify/verify.py that record how they were called.

    pytest tests/test_entrypoint.py
"""
import json
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
CHECK_RELEASE = ROOT / 'container' / 'check_release.py'

STUB = """#!/usr/bin/env python3
import json, os, sys
keys = ('PROHAP_HOME', 'PROHAP_REFS', 'ENSEMBL_RELEASE', 'PYTHONHASHSEED', 'PYTHONNOUSERSITE', 'HOME', 'XDG_CACHE_HOME')
print('STUB ' + json.dumps({'prog': os.path.basename(sys.argv[0]), 'argv': sys.argv[1:], 'cwd': os.getcwd(),
                            'env': {k: os.environ.get(k) for k in keys}}))
"""


@pytest.fixture
def image(tmp_path):
    """A stand-in for /opt/prohap with the references of release 113, and a bin/ with a stub snakemake."""
    home = tmp_path / 'opt' / 'prohap'
    (home / 'container').mkdir(parents=True)
    (home / 'verify').mkdir()
    (home / 'refs' / 'ensembl.v113').mkdir(parents=True)
    shutil.copy(ROOT / 'entrypoint.sh', home / 'entrypoint.sh')
    shutil.copy(CHECK_RELEASE, home / 'container' / 'check_release.py')
    (home / 'VERSION').write_text('9.9.9\n')
    (home / 'refs' / 'ensembl.v113' / 'release.txt').write_text('113\n')
    (home / 'verify' / 'verify.py').write_text(STUB)
    bindir = tmp_path / 'bin'
    bindir.mkdir()
    (bindir / 'snakemake').write_text(STUB)
    (bindir / 'snakemake').chmod(0o755)
    # python3 must be this interpreter (it has PyYAML); the entrypoint puts /opt/conda/envs/ProHap/bin first
    (bindir / 'python3').symlink_to(sys.executable)
    return home


@pytest.fixture
def work(tmp_path):
    work = tmp_path / 'work'
    work.mkdir()
    (work / 'config.yaml').write_text('ensembl_release: 113\nchromosomes: ["22"]\n')
    return work


def run(image, work, *args, env=None):
    environ = {'PATH': str(image.parent.parent / 'bin') + os.pathsep + '/usr/bin:/bin', 'HOME': str(work)}
    environ.update(env or {})
    return subprocess.run(['bash', str(image / 'entrypoint.sh')] + list(args), cwd=work, env=environ,
                          capture_output=True, text=True)


def stub_call(proc):
    lines = [line for line in proc.stdout.splitlines() if line.startswith('STUB ')]
    assert len(lines) == 1, proc.stdout + proc.stderr
    return json.loads(lines[0][len('STUB '):])


def test_no_arguments_is_a_dry_run(image, work):
    proc = run(image, work)
    assert proc.returncode == 0, proc.stderr
    call = stub_call(proc)
    assert call['prog'] == 'snakemake'
    assert call['argv'] == ['--snakefile', str(image / 'Snakefile'), '--directory', str(work),
                            '--rerun-triggers', 'mtime', '--dry-run']
    assert 'prohap-fast:9.9.9-ensembl.v113' in proc.stderr


@pytest.mark.parametrize('prefix', [[], ['snakemake']])
def test_arguments_go_to_snakemake(image, work, prefix):
    proc = run(image, work, *prefix, '--cores', '16', '-p')
    assert proc.returncode == 0, proc.stderr
    call = stub_call(proc)
    assert call['argv'][-3:] == ['--cores', '16', '-p']
    assert '--dry-run' not in call['argv']


def test_environment(image, work, tmp_path):
    """The references are the image's whatever the caller exported; the hash seed is fixed; no user site."""
    proc = run(image, work, '--cores', '1', env={'PROHAP_REFS': str(tmp_path / 'elsewhere'),
                                                  'ENSEMBL_RELEASE': '999', 'PYTHONHASHSEED': 'random'})
    assert proc.returncode == 0, proc.stderr
    env = stub_call(proc)['env']
    assert env['PROHAP_HOME'] == str(image)
    assert env['PROHAP_REFS'] == str(image / 'refs')
    assert env['ENSEMBL_RELEASE'] == '113'
    assert env['PYTHONHASHSEED'] == '0'
    assert env['PYTHONNOUSERSITE'] == '1'
    assert env['HOME'] == str(work)  # writable, so kept
    assert env['XDG_CACHE_HOME'] == str(work / '.cache')


def test_unwritable_home_moves_into_the_work_directory(image, work):
    proc = run(image, work, '--cores', '1', env={'HOME': '/nonexistent'})
    assert proc.returncode == 0, proc.stderr
    assert stub_call(proc)['env']['HOME'] == str(work / '.home')
    assert (work / '.home').is_dir()


def test_verify(image, work):
    proc = run(image, work, 'verify', '--cores', '2', '--transcripts', '20')
    assert proc.returncode == 0, proc.stderr
    call = stub_call(proc)
    assert call['prog'] == 'verify.py'
    assert call['argv'] == ['--cores', '2', '--transcripts', '20']
    assert call['cwd'] == str(work)
    assert call['env']['PYTHONHASHSEED'] == '0'


def test_prohap_work(image, work, tmp_path):
    other = tmp_path / 'other'
    other.mkdir()
    proc = run(image, other, 'verify', env={'PROHAP_WORK': str(work)})
    assert proc.returncode == 0, proc.stderr
    assert stub_call(proc)['cwd'] == str(work)


def test_shell(image, work):
    proc = run(image, work, 'shell', '-c', 'echo "$PYTHONNOUSERSITE $PROHAP_REFS"')
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == '1 %s\n' % (image / 'refs')


def test_help(image, work):
    proc = run(image, work, 'help')
    assert proc.returncode == 0
    assert 'verify [ARGS]' in proc.stdout and 'set -euo' not in proc.stdout
    assert proc.stdout.startswith('prohap-fast container entrypoint.')


def test_missing_config_names_the_template(image, work):
    (work / 'config.yaml').unlink()
    proc = run(image, work, '--cores', '1')
    assert proc.returncode == 1
    assert 'no config.yaml in the work directory' in proc.stderr
    assert str(image / 'config_tre_template.yaml') in proc.stderr
    assert 'STUB' not in proc.stdout


def test_other_release_stops_the_run(image, work):
    (work / 'config.yaml').write_text('ensembl_release: 114\n')
    for args in ([], ['verify']):
        proc = run(image, work, *args)
        assert proc.returncode == 1
        assert 'asks for Ensembl release 114' in proc.stderr
        assert 'carries the references of release 113' in proc.stderr
        assert 'prohap-fast:9.9.9-ensembl.v113' in proc.stderr
        assert 'STUB' not in proc.stdout


def test_unwritable_work_directory(image, work):
    if os.geteuid() == 0:
        pytest.skip('root can write anywhere')
    work.chmod(stat.S_IRUSR | stat.S_IXUSR)
    try:
        proc = run(image, work, '--cores', '1')
    finally:
        work.chmod(0o755)
    assert proc.returncode == 1
    assert 'is not writable' in proc.stderr


def test_image_without_references(image, work):
    shutil.rmtree(image / 'refs' / 'ensembl.v113')
    proc = run(image, work, '--cores', '1')
    assert proc.returncode == 1
    assert 'exactly one Ensembl release' in proc.stderr


# ------------------------------------ check_release.py ------------------------------------

def check_release(tmp_path, text, *args):
    config = tmp_path / 'config.yaml'
    config.write_text(text)
    return subprocess.run([sys.executable, str(CHECK_RELEASE), str(config)] + list(args), capture_output=True,
                          text=True)


@pytest.mark.parametrize('value', ['113', '"113"', "'113'"])
def test_check_release_accepts_the_image_release(tmp_path, value):
    proc = check_release(tmp_path, 'ensembl_release: %s\n' % value, '113')
    assert proc.returncode == 0, proc.stderr
    assert proc.stderr == ''


def test_check_release_refuses_another_release(tmp_path):
    proc = check_release(tmp_path, 'ensembl_release: 115\n', '113', 'prohap-fast:1.0.0-ensembl.v113')
    assert proc.returncode == 1
    assert 'asks for Ensembl release 115' in proc.stderr
    assert 'prohap-fast:1.0.0-ensembl.v113 carries the references of release 113' in proc.stderr
    assert 'ensembl_release: 113' in proc.stderr


def test_check_release_needs_the_key(tmp_path):
    proc = check_release(tmp_path, 'chromosomes: ["22"]\n', '113')
    assert proc.returncode == 1
    assert "has no 'ensembl_release'" in proc.stderr


@pytest.mark.parametrize('text', ['ensembl_release: [113\n', '- 113\n'])
def test_check_release_unreadable_config(tmp_path, text):
    assert check_release(tmp_path, text, '113').returncode == 2
