"""src/modules/fallback.py and src/modules/legacy_import.py:

    pytest tests/test_fast_helpers.py
"""
import os
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
FALLBACK_MESSAGE = 'ProHap-fast: falling back to the reference code'


def run_python(code, cwd, env=None):
    return subprocess.run([sys.executable, '-c', textwrap.dedent(code)], cwd=cwd, capture_output=True, text=True,
                          env=dict(os.environ, **(env or {})))


# ---------------------------------------------------------------------------------------------------------------
# legacy_import


def test_legacy_modules_come_from_legacy_src():
    result = run_python('''
        import sys
        sys.path.insert(0, 'src')
        import modules.common as own_common
        from modules.legacy_import import legacy
        import modules.common
        print(own_common.__file__)
        print(modules.common is own_common)
        for name in ('common', 'coordinates_toolbox', 'vcf_reader', 'get_haplotypes', 'process_haplotypes'):
            print(getattr(legacy, name).__file__)
        # the reference modules import each other, not src's
        print(legacy.process_haplotypes.get_rna_position is legacy.coordinates_toolbox.get_rna_position)
        print(legacy.vcf_reader.KeyWrapper is legacy.common.KeyWrapper)
        print(sorted(k for k in sys.modules if k == 'modules' or k.startswith('modules.')))
    ''', ROOT)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    assert Path(lines[0]).resolve() == ROOT / 'src' / 'modules' / 'common.py'
    assert lines[1] == 'True'
    assert [Path(p).resolve().parent for p in lines[2:7]] == [ROOT / 'legacy' / 'src' / 'modules'] * 5
    assert lines[7:9] == ['True', 'True']
    assert lines[9] == "['modules', 'modules.common', 'modules.legacy_import']"


def test_legacy_import_twice_is_the_same_object():
    result = run_python('''
        import sys
        sys.path.insert(0, 'src')
        from modules.legacy_import import legacy as a
        import importlib
        b = importlib.import_module('modules.legacy_import').legacy
        print(a is b, a.common is sys.modules['prohap_legacy.common'])
    ''', ROOT)
    assert result.returncode == 0, result.stderr
    assert result.stdout.split() == ['True', 'True']


# ---------------------------------------------------------------------------------------------------------------
# fallback


@pytest.fixture
def tree(tmp_path):
    """A stand-in repository: src/demo.py falls back to legacy/src/demo.py, which prints what it was given."""
    (tmp_path / 'src' / 'modules').mkdir(parents=True)
    (tmp_path / 'legacy' / 'src').mkdir(parents=True)
    (tmp_path / 'tools').mkdir()
    shutil.copy(ROOT / 'src' / 'modules' / 'fallback.py', tmp_path / 'src' / 'modules')
    shutil.copy(ROOT / 'tools' / 'run_fork.py', tmp_path / 'tools')
    (tmp_path / 'src' / 'demo.py').write_text(textwrap.dedent('''
        import sys
        from modules.fallback import to_legacy
        print('fast started', flush=True)
        to_legacy(__file__, 'demo reason')
        print('not reached')
    '''))
    (tmp_path / 'legacy' / 'src' / 'demo.py').write_text(textwrap.dedent('''
        import multiprocessing, os, sys
        print('legacy', sys.argv[1:], multiprocessing.get_start_method(), os.environ.get('PYTHONHASHSEED'))
    '''))
    return tmp_path


def test_fallback_runs_the_legacy_script_with_the_same_arguments(tree):
    result = subprocess.run([sys.executable, str(tree / 'src' / 'demo.py'), '-i', 'a b', '-o', 'out'],
                            capture_output=True, text=True, env=dict(os.environ, PYTHONHASHSEED='0'))
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ['fast started', "legacy ['-i', 'a b', '-o', 'out'] fork 0"]
    assert result.stderr.strip() == FALLBACK_MESSAGE + ' (demo.py): demo reason'


def test_fallback_through_the_fork_wrapper(tree):
    """As the Snakefile runs scripts on macOS: tools/run_fork.py SCRIPT ARGS."""
    result = subprocess.run([sys.executable, str(tree / 'tools' / 'run_fork.py'), str(tree / 'src' / 'demo.py'),
                             'x'], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines()[-1].startswith("legacy ['x'] fork")


def test_fallback_error_mode(tree):
    result = subprocess.run([sys.executable, str(tree / 'src' / 'demo.py')], capture_output=True, text=True,
                            env=dict(os.environ, PROHAP_FALLBACK='error'))
    assert result.returncode == 3
    assert result.stdout == 'fast started\n'
    assert FALLBACK_MESSAGE in result.stderr


def test_every_src_script_has_a_legacy_counterpart():
    scripts = sorted(p.name for p in (ROOT / 'src').glob('*.py'))
    assert scripts == sorted(p.name for p in (ROOT / 'legacy' / 'src').glob('*.py'))
