"""The filter fused into prohap.py (F3, IMPLEMENTATION_NOTES section 8.2).

Without fusion the workflow runs vcf_filter_fix.py, which writes the filtered VCF gzipped, and prohap.py, which
decompresses and reads it again. With `prohap.py -filter_glob GLOB -filter_af AF -filter_af_field FIELD` (instead of
`-i FILE`), prohap.py starts vcf_filter_fix.py with `-o -` and reads its output through a pipe: the same lines, so the
same line numbers for `.` IDs, without writing or compressing them. The filter runs in its own process, in parallel
with prohap.py's reader.

Input the fast filter does not handle stops it with exit status 3 (it cannot hand a pipe to the reference); input
prohap.py does not handle raises Unsupported as before. Either way `to_legacy_fused` runs the two steps separately:
vcf_filter_fix.py (which itself falls back if it must) into a file in -tmp_dir, then the reference prohap.py on it.
"""
import os
import subprocess
import sys
from pathlib import Path

from modules.fallback import ERROR_STATUS, MESSAGE, Unsupported, legacy_command

FILTER = Path(__file__).resolve().parent.parent / 'vcf_filter_fix.py'
FUSED_OPTIONS = ('-filter_glob', '-filter_af', '-filter_af_field')


def filter_command(args, output):
    return [sys.executable, str(FILTER), '-i', args.filter_glob, '-chr', args.chromosome, '-af', str(args.filter_af),
            '-af_field', args.filter_af_field, '-o', output]


class FilterReader:
    """The lines of the filtered VCF as bytes, read from vcf_filter_fix.py -o - (the interface of sweep.VCFReader)."""

    def __init__(self, args):
        env = dict(os.environ, PROHAP_FALLBACK='error')
        self.process = subprocess.Popen(filter_command(args, '-'), stdout=subprocess.PIPE, bufsize=1 << 20, env=env)
        self.handle = self.process.stdout
        self.path = args.filter_glob

    def readline(self):
        line = self.handle.readline()
        if b'\r' in line or not line.isascii():
            raise Unsupported('carriage return or non-ASCII byte in the filtered VCF')
        return line

    def close(self):
        self.handle.close()
        if self.process.poll() is None:
            self.process.kill()
        self.process.wait()

    def finish(self):
        """Close after reading everything; the filter's exit status says whether everything was filtered."""
        self.handle.close()
        code = self.process.wait()
        if code == ERROR_STATUS:
            raise Unsupported('the filter met input it does not handle (see above)')
        if code != 0:
            raise RuntimeError('vcf_filter_fix.py failed (exit status %d) on %s' % (code, self.path))


def unfused_argv(argv, filtered):
    """argv with the -filter_* options replaced by -i FILTERED."""
    out = []
    skip = False
    for arg in argv:
        if skip:
            skip = False
            continue
        if arg in FUSED_OPTIONS:
            skip = True
            continue
        if any(arg.startswith(o + '=') for o in FUSED_OPTIONS):
            continue
        out.append(arg)
    return ['-i', filtered] + out


def to_legacy_fused(script, args, reason, argv=None):
    """Run the filter into a file, then the reference prohap.py on it (with argv, default this process's arguments,
    the -filter_* options replaced by -i), and exit with its status. Does not return."""
    sys.stdout.flush()
    print('%s (%s): %s' % (MESSAGE, Path(script).name, reason), file=sys.stderr, flush=True)
    if os.environ.get('PROHAP_FALLBACK') == 'error':
        print('PROHAP_FALLBACK=error: stopping instead', file=sys.stderr, flush=True)
        os._exit(ERROR_STATUS)
    os.makedirs(args.tmp_dir, exist_ok=True)
    filtered = os.path.join(args.tmp_dir, 'chr%s_phased_filtered.vcf.gz' % args.chromosome)
    code = subprocess.call(filter_command(args, filtered))
    if code == 0:
        code = subprocess.call(legacy_command(script, unfused_argv(sys.argv[1:] if argv is None else argv,
                                                                   filtered)))
    if os.path.exists(filtered):
        os.remove(filtered)
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(code)
