"""Hand a run over to the reference code.

ProHap-fast never guesses at input it does not support. When one of its scripts meets such input, it replaces its own
process with the script of the same name in legacy/src, run with the same arguments by the same interpreter. The
reference then writes every output itself, from the start, overwriting anything the fast script may have written.

    from modules.fallback import to_legacy
    to_legacy(__file__, 'unphased genotype at line 1234')

The reason goes to stderr, after MESSAGE, so a test can tell that a run fell back. With PROHAP_FALLBACK=error in the
environment the script stops with exit status 3 instead (tests use it to make sure the fast path ran).

Before calling to_legacy, a script must stop its worker processes (Pool.terminate()): forked workers would outlive the
exec and wait for work forever.
"""
import multiprocessing
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LEGACY_SRC = ROOT / 'legacy' / 'src'
RUN_FORK = ROOT / 'tools' / 'run_fork.py'
MESSAGE = 'ProHap-fast: falling back to the reference code'
ERROR_STATUS = 3


class Unsupported(Exception):
    """Raised by fast code for input it does not handle; the script catches it and calls to_legacy."""


def legacy_command(script, argv):
    """The command that runs legacy/src/<name of script> with argv, as this interpreter would run it."""
    target = LEGACY_SRC / Path(script).name
    if not target.is_file():
        raise FileNotFoundError('no reference script ' + str(target))
    command = [sys.executable] + subprocess._args_from_interpreter_flags()
    # The reference hands closures to multiprocessing.Pool, which needs workers started by fork. Where fork is not
    # the default (macOS), run it through the wrapper that selects it, as the Snakefile does.
    if multiprocessing.get_all_start_methods()[0] != 'fork':
        command.append(str(RUN_FORK))
    return command + [str(target)] + list(argv)


def cannot_fall_back(script, reason):
    """Stop with ERROR_STATUS where to_legacy would hand the run over, because the reference cannot write the output
    asked for (prohap.py -fix_bugs 1: the reference has none of the fixes). Does not return."""
    sys.stdout.flush()
    print('ProHap-fast: this input needs the reference code (%s): %s' % (Path(script).name, reason), file=sys.stderr)
    print('The reference code has none of the fix_bugs fixes, so it cannot write this output. Set fix_bugs: false '
          '(prohap.py -fix_bugs 0) to run it with the reference code, or change the input.', file=sys.stderr,
          flush=True)
    os._exit(ERROR_STATUS)


def to_legacy(script, reason, argv=None):
    """Replace this process with the reference script of the same name and the same arguments. Does not return."""
    command = legacy_command(script, sys.argv[1:] if argv is None else argv)
    sys.stdout.flush()
    print('%s (%s): %s' % (MESSAGE, Path(script).name, reason), file=sys.stderr, flush=True)
    if os.environ.get('PROHAP_FALLBACK') == 'error':
        print('PROHAP_FALLBACK=error: stopping instead', file=sys.stderr, flush=True)
        os._exit(ERROR_STATUS)
    os.execv(command[0], command)
