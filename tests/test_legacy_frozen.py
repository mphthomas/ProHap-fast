"""legacy/src must be upstream ProHap `cafcc51` src/ plus legacy/legacy-fixes.patch, byte for byte.

legacy/src is the frozen reference that ProHap-fast's output is checked against (see legacy/README.md). This test
rebuilds it from git history and compares every file, symlink and executable bit. It needs the full history: in
GitHub Actions, use actions/checkout with `fetch-depth: 0`.

    pytest tests/test_legacy_frozen.py
    python tests/test_legacy_frozen.py                 # same check without pytest
    python tests/test_legacy_frozen.py --write-patch   # after a deliberate output change to legacy/src
"""
import io
import os
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LEGACY = ROOT / 'legacy'
LEGACY_SRC = LEGACY / 'src'
PATCH = LEGACY / 'legacy-fixes.patch'
UPSTREAM_COMMIT = 'cafcc515ba3ba6b1b6121db6ea3ac612ca053124'  # full SHA: a ref named like the short hash cannot stand in
SKIP_DIRS = {'__pycache__'}  # untracked bytecode written when the legacy arm runs


def _git_env(top):
    """Environment in which git sees only the repository at `top` (or none), whatever the caller exported."""
    env = {k: v for k, v in os.environ.items() if not k.startswith('GIT_')}
    env['GIT_CEILING_DIRECTORIES'] = str(Path(top).resolve().parent)
    return env


def _git(args, cwd, env, **kwargs):
    return subprocess.run(['git'] + args, cwd=cwd, env=env, capture_output=True, **kwargs)


def extract_upstream(dest):
    """Extract `git archive cafcc51 src` into dest."""
    env = _git_env(ROOT)
    archive = _git(['--no-replace-objects', 'archive', UPSTREAM_COMMIT, 'src'], ROOT, env)
    if archive.returncode != 0:
        hint = ''
        if _git(['rev-parse', '--is-shallow-repository'], ROOT, env, text=True).stdout.strip() == 'true':
            hint = '\nThis is a shallow clone; the check needs the full history (actions/checkout: fetch-depth: 0).'
        raise RuntimeError(f'git archive {UPSTREAM_COMMIT[:7]} failed: {archive.stderr.decode().strip()}{hint}')
    with tarfile.open(fileobj=io.BytesIO(archive.stdout)) as tar:
        tar.extractall(dest, **({'filter': 'data'} if hasattr(tarfile, 'data_filter') else {}))


def build_reference(dest):
    """cafcc51 src/ with legacy-fixes.patch applied, in dest; returns dest/src."""
    extract_upstream(dest)
    # dest is outside any repository, so git apply works as a plain patch tool there
    applied = _git(['apply', '--whitespace=nowarn', str(PATCH)], dest, _git_env(dest))
    if applied.returncode != 0:
        raise RuntimeError(f'legacy-fixes.patch does not apply to {UPSTREAM_COMMIT[:7]}: '
                           + applied.stderr.decode().strip())
    return Path(dest) / 'src'


def tree_files(top):
    """{path relative to top: (content, executable)}. Symlinks, to files or directories, are recorded as links and
    never followed; __pycache__ directories are skipped."""
    files = {}
    for dirpath, dirnames, filenames in os.walk(top):
        for d in list(dirnames):
            if (Path(dirpath) / d).is_symlink():  # os.walk lists these as directories but never enters them
                dirnames.remove(d)
                filenames.append(d)
            elif d in SKIP_DIRS:
                dirnames.remove(d)
        for name in filenames:
            path = Path(dirpath) / name
            mode = path.lstat().st_mode
            content = (b'symlink -> ' + os.readlink(path).encode()) if stat.S_ISLNK(mode) else path.read_bytes()
            files[path.relative_to(top).as_posix()] = (content, bool(mode & stat.S_IXUSR))
    return files


def _first_differing_line(a, b):
    i = next((k for k, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))
    return a[:i].count(b'\n') + 1


def differences(expected_top):
    """Human-readable differences between the rebuilt reference and legacy/src; empty if identical."""
    out = [f'{p.relative_to(ROOT)} must be a real directory, not a symlink'
           for p in (LEGACY, LEGACY_SRC) if not stat.S_ISDIR(p.lstat().st_mode)]
    if out:
        return out
    tracked = _git(['ls-files', '-z', '--', 'legacy/src'], ROOT, _git_env(ROOT)).stdout.decode().split('\0')
    out += [f'bytecode tracked by git: {p}' for p in tracked if SKIP_DIRS & set(p.split('/'))]
    expected, actual = tree_files(expected_top), tree_files(LEGACY_SRC)
    out += [f'missing from legacy/src: {p}' for p in sorted(expected.keys() - actual.keys())]
    out += [f'not in {UPSTREAM_COMMIT[:7]} + patch: legacy/src/{p}' for p in sorted(actual.keys() - expected.keys())]
    for p in sorted(expected.keys() & actual.keys()):
        (e_bytes, e_exec), (a_bytes, a_exec) = expected[p], actual[p]
        if e_bytes != a_bytes:
            out.append(f'content differs: legacy/src/{p}, first at line {_first_differing_line(e_bytes, a_bytes)}')
        elif e_exec != a_exec:
            out.append(f'executable bit differs: legacy/src/{p}')
    return out


def write_patch():
    """Regenerate legacy-fixes.patch as the full diff from cafcc51 src/ to the current legacy/src, independent of the
    caller's git configuration."""
    with tempfile.TemporaryDirectory() as tmp:
        extract_upstream(tmp)
        env = _git_env(tmp)
        env.update(GIT_CONFIG_NOSYSTEM='1', GIT_CONFIG_GLOBAL=os.devnull, HOME=tmp, XDG_CONFIG_HOME=tmp)
        git = ['-c', 'user.name=legacy', '-c', 'user.email=legacy@invalid', '-c', 'core.autocrlf=false',
               '-c', 'core.abbrev=7', '-c', 'diff.algorithm=myers']
        for args in (['init', '-q'], ['add', '-A', '-f'], ['commit', '-q', '-m', 'cafcc51 src']):
            _git(git + args, tmp, env, check=True)
        shutil.rmtree(Path(tmp) / 'src')
        shutil.copytree(LEGACY_SRC, Path(tmp) / 'src', symlinks=True, ignore=shutil.ignore_patterns(*SKIP_DIRS))
        _git(git + ['add', '-A', '-f', '--', 'src'], tmp, env, check=True)
        diff = _git(git + ['diff', '--cached', '--no-color', '--no-ext-diff', '--no-renames', '--indent-heuristic',
                           '-U3', '--', 'src'], tmp, env, check=True).stdout
    PATCH.write_bytes(diff)


def test_legacy_src_is_cafcc51_plus_patch(tmp_path):
    problems = differences(build_reference(tmp_path))
    assert not problems, 'legacy/src is not cafcc51 + legacy/legacy-fixes.patch:\n  ' + '\n  '.join(problems)


if __name__ == '__main__':
    if '--write-patch' in sys.argv[1:]:
        write_patch()
        print(f'wrote {PATCH.relative_to(ROOT)}')
    with tempfile.TemporaryDirectory() as tmp:
        problems = differences(build_reference(tmp))
    print('\n'.join(problems) if problems else 'OK: legacy/src == cafcc51 + legacy/legacy-fixes.patch')
    sys.exit(1 if problems else 0)
