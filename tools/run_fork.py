"""Run one of ProHap's scripts with the multiprocessing start method set to 'fork'.

ProHap hands closures to multiprocessing.Pool, which works only when workers are forked, the default on Linux.
macOS defaults to 'spawn', so there the Snakefile runs every script through this wrapper.

Usage: python3 tools/run_fork.py SCRIPT [ARGS...]    (behaves like: python3 SCRIPT [ARGS...])
"""
import multiprocessing
import os
import runpy
import sys

if len(sys.argv) < 2:
    sys.exit(__doc__)
multiprocessing.set_start_method('fork')
script = sys.argv[1]
sys.argv = sys.argv[1:]
# like `python3 SCRIPT`: the script's own directory comes first, so `from modules...` resolves in its tree
sys.path[0] = os.path.dirname(os.path.realpath(script))
runpy.run_path(script, run_name='__main__')
