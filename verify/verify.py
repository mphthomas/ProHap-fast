#!/usr/bin/env python3
"""verify: check ProHap-fast against the frozen reference implementation on your own data.

Run it from the directory that holds your config.yaml, the directory a normal run uses:

    python3 <ProHap-fast>/verify/verify.py --cores 8 [--chromosomes 21,22]
        [--samples 2000 [--seed 0] | --sample-list ids.txt] [--transcripts 50] [--keep] [--dry-run]

verify runs the workflow twice, in two work directories, verify_runs/<timestamp>/legacy/ and
verify_runs/<timestamp>/fast/, each with a copy of config.yaml that differs only in `implementation` (legacy:
legacy/src, the frozen reference; fast: src/). Both read the same inputs and references; relative input paths are
resolved against the current directory, and every output path is moved into the arm's work directory, so nothing
outside verify_runs/ is written. Then it compares every output, decompressed, byte for byte, and prints a table
(file, both sha256s, IDENTICAL or DIFFERENT, first differing line) and both run times. The exit status is 0 only if
every output is identical, 1 if any differs, 2 if verify or the legacy run failed.

On success the work directories are deleted and only report.tsv is kept; --keep keeps everything. --dry-run
prepares the run directory and dry-runs both arms.

Subsetting, applied to both arms alike:
  --chromosomes 21,22    only these chromosomes of the config's list
  --samples N            N samples drawn at random (seeded by --seed, default 0), kept in VCF header order;
  --sample-list FILE     or the samples named in FILE, one per line. The VCF chunks are cut to those samples,
                         changing nothing else (line numbers name the records whose ID is "."), and
                         sample_metadata_file to the matching rows.
  --transcripts N        only the first N transcripts of each chromosome in the transcript list. Keeps every sample,
                         which matters: the reference reads its tables in pandas chunks whose size depends on the
                         number of samples.

Each arm needs its own work directory: switching `implementation` in a directory that already holds outputs reruns
nothing, because Snakemake sees the same rules with the same code and the outputs up to date.

A fast script that meets input it does not handle hands it to the reference script of the same name, which then
writes that step's outputs; verify lists such steps after the fast arm has run. Their outputs are compared all the
same, but for them the comparison checks the reference against itself.
"""
import argparse
import concurrent.futures
import datetime
import glob
import os
import random
import re
import resource
import shlex
import shutil
import subprocess
import sys
import time

import yaml

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import compare  # noqa: E402
import preflight  # noqa: E402

ROOT = os.path.dirname(HERE)
SNAKEFILE = os.path.join(ROOT, 'Snakefile')
ARMS = ('legacy', 'fast')

# Paths the workflow writes. Relative ones stay as they are, inside the arm's work directory.
OUTPUT_KEYS = ('haplo_table_file', 'haplo_fasta_file', 'haplo_cdna_file', 'final_fasta_file',
               'var_table_file', 'var_fasta_file', 'var_cdna_file')
# Paths the workflow reads, relative to the work directory. phased_local_path, custom_transcript_list, refs_dir
# and variant_vcf are handled separately.
INPUT_KEYS = ('contaminants_fasta', 'sample_metadata_file', 'haplo_added_table', 'haplo_added_fasta')


class VerifyError(Exception):
    pass


def say(*args):
    print('[verify]', *args, flush=True)


def parse_args(argv):
    p = argparse.ArgumentParser(prog='verify', description=__doc__.split('\n\n')[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog='\n\n'.join(__doc__.split('\n\n')[1:]))
    p.add_argument('--cores', type=int, default=1, help='cores for each Snakemake run (default 1)')
    p.add_argument('--chromosomes', help='comma-separated subset of the config\'s chromosomes, e.g. 21,22')
    g = p.add_mutually_exclusive_group()
    g.add_argument('--samples', type=int, metavar='N', help='use N samples drawn at random')
    g.add_argument('--sample-list', metavar='FILE', help='use the samples named in FILE, one per line')
    p.add_argument('--seed', type=int, default=0, help='random seed for --samples (default 0)')
    p.add_argument('--transcripts', type=int, metavar='N', help='use the first N transcripts of each chromosome')
    p.add_argument('--keep', action='store_true', help='keep the run directory when every output is identical')
    p.add_argument('--dry-run', action='store_true', help='prepare the run directory and dry-run both arms')
    # development: run the fast arm with the scripts in DIR instead of src/
    p.add_argument('--fast-src', metavar='DIR', help=argparse.SUPPRESS)
    args = p.parse_args(argv)
    if args.cores < 1:
        p.error('--cores must be at least 1')
    if args.samples is not None and args.samples < 1:
        p.error('--samples must be at least 1')
    if args.transcripts is not None and args.transcripts < 1:
        p.error('--transcripts must be at least 1')
    if args.fast_src is not None:
        if not os.path.isfile(os.path.join(args.fast_src, 'prohap.py')):
            p.error('--fast-src %s has no prohap.py' % args.fast_src)
        args.fast_src = os.path.abspath(args.fast_src)
    return args


# ------------------------------------ config ------------------------------------

def load_config(path):
    try:
        with open(path) as fh:
            cfg = yaml.safe_load(fh)
    except OSError as e:
        raise VerifyError('cannot read %s: %s. Run verify from the directory that holds your config.yaml.'
                          % (path, e.strerror))
    if not isinstance(cfg, dict):
        raise VerifyError(path + ' is not a YAML mapping')
    missing = [k for k in ('ensembl_release', 'chromosomes', 'phased_vcf_file_name', 'sample_metadata_file')
               if k not in cfg]
    if missing:
        raise VerifyError('%s has no %s' % (path, ', '.join(repr(k) for k in missing)))
    if not cfg.get('use_ProHap'):
        raise VerifyError('use_ProHap is off in %s, so there is nothing to compare: the ProVar rules run the legacy '
                          'code in both arms' % path)
    return cfg


def refs_root(cfg):
    """The references directory exactly as the Snakefile resolves it, from the current directory."""
    return os.path.abspath(os.path.expanduser(str(cfg.get('refs_dir') or os.environ.get('PROHAP_REFS')
                                                  or '/opt/prohap/refs')))


def phased_pattern(cfg):
    """The VCF glob of the filter rule, before {chr} is filled in, relative to the work directory."""
    local = str(cfg.get('phased_local_path') or '')
    return (local if local else 'data/vcf/phased/') + str(cfg['phased_vcf_file_name'])


def chunks(pattern, chrom, restrict=False):
    """The VCF chunks of one chromosome, in the order the legacy filter reads them (restrict_to_exons: index files
    left out, as tools/restrict_to_exons.py does)."""
    paths = sorted(glob.glob(pattern.replace('{chr}', str(chrom))))
    return [p for p in paths if not p.endswith(('.tbi', '.csi'))] if restrict else paths


def transcript_source(cfg, refs):
    """Where the Snakefile takes the transcript list from: ('file', path) or ('rule', work-dir path built by a rule).
    Mirrors TRANSCRIPT_LIST in the Snakefile."""
    release = str(cfg['ensembl_release'])
    workdir_list = 'data/transcripts_reference_' + release + '.csv'
    mane_list = os.path.join(refs, 'ensembl.v' + release, 'transcripts_reference_' + release + '.csv')
    reference = ('file', mane_list) if cfg.get('only_MANE_select') else ('rule', workdir_list)
    custom = cfg.get('custom_transcript_list')
    if custom is None:
        return reference
    if len(custom) == 0:
        return ('rule', 'data/included_transcripts.csv')
    if os.path.normpath(custom) == workdir_list:
        return reference
    return ('file', custom)


class Run:
    """One verify run: the run directory and the config of each arm."""

    def __init__(self, cfg, args, cwd):
        self.cfg = cfg
        self.args = args
        self.cwd = cwd
        stamp = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        base = os.path.join(cwd, 'verify_runs', stamp)
        run_dir, n = base, 1
        while os.path.exists(run_dir):
            n += 1
            run_dir = '%s-%d' % (base, n)
        self.dir = run_dir
        self.input_dir = os.path.join(run_dir, 'input')
        # the arms' work directories (and prep/) sit three levels below cwd: verify_runs/<stamp>/<arm>
        self.up = os.path.join('..', '..', '..')
        self.refs = refs_root(cfg)
        self.pattern = phased_pattern(cfg)
        self.chromosomes = self.select_chromosomes()
        self.chunks = {c: chunks(self.pattern, c, bool(cfg.get('restrict_to_exons'))) for c in self.chromosomes}
        self.overrides = {}

    def arm_dir(self, arm):
        return os.path.join(self.dir, arm)

    def from_arm(self, path):
        """An input path given relative to cwd, as seen from an arm's work directory."""
        return path if os.path.isabs(path) else os.path.join(self.up, path)

    def select_chromosomes(self):
        listed = self.cfg.get('chromosomes') or []
        if not self.args.chromosomes:
            return list(listed)
        wanted = [c.strip() for c in self.args.chromosomes.split(',') if c.strip()]
        names = [str(c) for c in listed]
        unknown = [c for c in wanted if c not in names]
        if unknown:
            raise VerifyError('--chromosomes %s: not in the config\'s chromosomes (%s)'
                              % (','.join(unknown), ', '.join(names)))
        # the config's order: the merge concatenates the chromosomes in that order
        return [c for c in listed if str(c) in wanted]

    def arm_config(self, arm):
        """The config of one arm: the user's, with implementation set, inputs resolved and outputs moved."""
        cfg = dict(self.cfg)
        cfg['implementation'] = arm
        cfg['refs_dir'] = self.refs
        cfg['chromosomes'] = list(self.chromosomes)
        for key in INPUT_KEYS:
            if cfg.get(key):
                cfg[key] = self.from_arm(str(cfg[key]))
        local = str(cfg.get('phased_local_path') or '') or 'data/vcf/phased/'
        cfg['phased_local_path'] = self.from_arm(local)
        custom = cfg.get('custom_transcript_list')
        if custom and transcript_source(self.cfg, self.refs) == ('file', custom):
            cfg['custom_transcript_list'] = self.from_arm(str(custom))
        if isinstance(cfg.get('variant_vcf'), dict):
            cfg['variant_vcf'] = {name: dict(entry, file=self.from_arm(str(entry['file'])))
                                  if isinstance(entry, dict) and entry.get('file') else entry
                                  for name, entry in cfg['variant_vcf'].items()}
        for key in OUTPUT_KEYS:
            if cfg.get(key):
                cfg[key] = self.arm_output(str(cfg[key]))
        cfg.update(self.overrides)
        return cfg

    def arm_output(self, path):
        """An output path inside the arm's work directory. Relative paths stay as they are; absolute ones, and
        relative ones that leave the work directory, move under outputs/ with their full path."""
        norm = os.path.normpath(path)
        if not os.path.isabs(norm) and norm != '..' and not norm.startswith('..' + os.sep):
            return path
        return os.path.join('outputs', os.path.abspath(os.path.join(self.cwd, path)).lstrip(os.sep))

    def write_configs(self):
        for arm in ARMS:
            os.makedirs(self.arm_dir(arm), exist_ok=True)
            write_yaml(os.path.join(self.arm_dir(arm), 'config.yaml'), self.arm_config(arm))


def write_yaml(path, data):
    with open(path, 'w') as fh:
        yaml.safe_dump(data, fh, sort_keys=False, default_flow_style=None, width=120)


# ------------------------------------ subsetting ------------------------------------

def read_sample_list(path):
    try:
        with open(path) as fh:
            names = [line.strip() for line in fh]
    except OSError as e:
        raise VerifyError('--sample-list %s: %s' % (path, e.strerror))
    names = [n for n in names if n and not n.startswith('#')]
    if not names:
        raise VerifyError('--sample-list %s names no samples' % path)
    return names


def select_samples(run):
    """The samples to keep, in the header order of the first chunk of the first chromosome; None to keep all."""
    if run.args.samples is None and run.args.sample_list is None:
        return None
    first = next((paths[0] for paths in run.chunks.values() if paths), None)
    if first is None:
        raise VerifyError('no VCF matches ' + run.pattern + ' for any chromosome, so there are no samples to pick from')
    header = preflight.read_samples(first)
    if not header:
        raise VerifyError(first + ' has no #CHROM header line')
    if run.args.sample_list is not None:
        wanted = read_sample_list(run.args.sample_list)
        unknown = sorted(set(wanted) - set(header))
        if unknown:
            raise VerifyError('--sample-list: %d samples are not in %s: %s'
                              % (len(unknown), first, ', '.join(unknown[:5]) + (', ...' if len(unknown) > 5 else '')))
        wanted = set(wanted)
        return [s for s in header if s in wanted]
    n = run.args.samples
    if n > len(header):
        raise VerifyError('--samples %d: %s has only %d samples' % (n, first, len(header)))
    if n == len(header):
        say('--samples %d is every sample: the VCFs are used as they are' % n)
        return None
    keep = sorted(random.Random(run.args.seed).sample(range(len(header)), n))
    return [header[i] for i in keep]


def subset_vcfs(run, samples):
    """Cut every chunk to `samples` with vcfcut.py into input/vcf/chr<c>/, several chunks at a time. The cut files
    differ from the originals only by the columns of the other samples, so every record keeps its line number. Chunk
    names are kept (prefixed with their position when a chromosome's chunks come from more than one directory), so
    the filter reads them in the same order and as gzip or text, as before."""
    restrict = bool(run.cfg.get('restrict_to_exons'))
    samples_file = os.path.join(run.input_dir, 'samples.txt')
    with open(samples_file, 'w') as fh:
        fh.write(''.join(s + '\n' for s in samples))
    jobs = []
    for chrom, paths in run.chunks.items():
        out_dir = os.path.join(run.input_dir, 'vcf', 'chr' + str(chrom))
        os.makedirs(out_dir, exist_ok=True)
        same_dir = len({os.path.dirname(p) for p in paths}) <= 1
        for i, path in enumerate(paths):
            name = os.path.basename(path) if same_dir else '%04d-%s' % (i, os.path.basename(path))
            dst = os.path.join(out_dir, name)
            if restrict and path.endswith('.bcf'):
                raise VerifyError('%s is BCF: --samples and --sample-list cut text VCF only. Convert the chunks to '
                                  'bgzip-compressed VCF (bcftools view -Oz), or verify with every sample' % path)
            if restrict and dst.endswith('.gz'):
                # restrict_to_exons reads indexed chunks: cut to text, then bgzip and index
                cmd = '%s %s %s %s %s && bgzip -f %s && tabix -f -p vcf %s' % tuple(shlex.quote(a) for a in (
                    sys.executable, os.path.join(HERE, 'vcfcut.py'), path, dst[:-3], samples_file, dst[:-3], dst))
                jobs.append((path, ['bash', '-c', 'set -euo pipefail; ' + cmd]))
            else:
                jobs.append((path, [sys.executable, os.path.join(HERE, 'vcfcut.py'), path, dst, samples_file]))

    def cut(job):
        return job[0], subprocess.run(job[1], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                                      env=child_env())
    with concurrent.futures.ThreadPoolExecutor(max_workers=run.args.cores) as pool:
        failed = [(path, proc.stderr.strip()) for path, proc in pool.map(cut, jobs) if proc.returncode != 0]
    if failed:
        raise VerifyError('could not cut %s to the selected samples:\n%s' % (failed[0][0], failed[0][1]))
    run.overrides['phased_local_path'] = os.path.join('..', 'input', 'vcf') + os.sep
    run.overrides['phased_vcf_file_name'] = 'chr{chr}' + os.sep + '*'


def read_metadata(path):
    """The metadata table as prohap.py reads it."""
    import pandas as pd
    return pd.read_csv(path, sep='\t')


def subset_metadata(run, samples):
    """sample_metadata_file cut to the rows of `samples`. Lines are copied as they are when the file is plain text
    whose lines map one to one onto pandas' rows; otherwise the rows are written back through pandas."""
    src = str(run.cfg['sample_metadata_file'])
    df = read_metadata(src)
    if 'Sample name' not in df.columns:
        raise VerifyError("sample_metadata_file %s has no 'Sample name' column" % src)
    keep = df['Sample name'].isin(set(samples)).tolist()
    dst = os.path.join(run.input_dir, os.path.basename(src))
    copied = False
    if not re.search(r'\.(gz|bz2|zip|xz|zst|tar)$', src):
        with open(src, 'rb') as fh:
            lines = fh.read().split(b'\n')
        trailing = lines and lines[-1] == b''
        if trailing:
            lines.pop()
        header, rows = lines[0], [line for line in lines[1:] if line.strip(b'\r') != b'']
        if len(rows) == len(keep) and b'"' not in b'\n'.join(lines):
            with open(dst, 'wb') as fh:
                fh.write(b'\n'.join([header] + [r for r, k in zip(rows, keep) if k]) + (b'\n' if trailing else b''))
            copied = True
    if not copied:
        df[keep].to_csv(dst, sep='\t', index=False)
    run.overrides['sample_metadata_file'] = os.path.join('..', 'input', os.path.basename(src))
    return sum(keep)


def subset_transcripts(run, n, prepare_list):
    """The first n transcripts of each chromosome of the transcript list, as prohap.py selects them, written to
    input/ and set as custom_transcript_list. prepare_list builds a list that a rule makes in the work directory."""
    import pandas as pd
    kind, path = transcript_source(run.cfg, run.refs)
    if kind == 'rule':
        path = prepare_list(path)
    try:
        df = pd.read_csv(path)
    except OSError as e:
        raise VerifyError('cannot read the transcript list %s: %s' % (path, e.strerror))
    chrom = df['chromosome'].apply(lambda x: str(x))
    keep = pd.Series(False, index=df.index)
    for c in run.chromosomes:
        rows = chrom[chrom == str(c)].index[:n]
        keep[rows] = True
    dst = os.path.join(run.input_dir, 'transcripts_first%d.csv' % n)
    df[keep].to_csv(dst, index=False)
    run.overrides['custom_transcript_list'] = os.path.join('..', 'input', os.path.basename(dst))
    return {str(c): int((chrom[keep] == str(c)).sum()) for c in run.chromosomes}


# ------------------------------------ outputs ------------------------------------

def strip_ext(path, n):
    return '.'.join(path.split('.')[:-n])


def expected_outputs(cfg, refs):
    """(path in the work directory, required) for every file the scripts write, in workflow order; computed from an
    arm's config as the Snakefile names them. Files that are not required may be missing from both arms."""
    out = []
    release = str(cfg['ensembl_release'])
    kind, path = transcript_source(cfg, refs)
    if kind == 'rule':
        out.append((path, False))
    out.append(('data/fasta/ensembl_reference_proteinDB_' + release + '_clean.fa.gz', True))
    out.append(('data/fasta/crap_tagged.fa.gz', False))     # only when contaminants_fasta is not ProHap's
    haplo = 'results/' + str(cfg['working_dir_name_haplo'])
    fasta_gz = '.gz' if str(cfg['haplo_fasta_file']).endswith('.gz') else ''
    for c in cfg['chromosomes']:
        if cfg.get('restrict_to_exons'):
            out.append(('data/vcf/phased/chr%s_exons.vcf.gz' % c, True))
        out.append(('data/vcf/phased/chr%s_phased_filtered.vcf.gz' % c, True))
        out.append(('%s/haplo_chr%s.tsv.gz' % (haplo, c), True))
        out.append(('%s/haplo_chr%s.fa%s' % (haplo, c, fasta_gz), True))
        out.append(('log/prohap_chr%s.log' % c, False))      # the removed-samples log; none for an empty VCF
    out.append((haplo + '/haplo_all.tsv.gz', True))
    table = str(cfg['haplo_table_file'])
    out.append((table, True))
    out.append(((strip_ext(table, 2) + '_sampleIDs.tsv.gz') if table.endswith('.gz')
                else (strip_ext(table, 1) + '_sampleIDs.tsv'), True))
    out.append((str(cfg['haplo_fasta_file']), True))
    if cfg.get('haplo_cdna_file'):
        out.append((str(cfg['haplo_cdna_file']), True))
    out.append(('results/haplo_all_clean.fa.gz', True))
    if cfg.get('use_ProVar'):
        out += [(str(cfg['var_table_file']), True), (str(cfg['var_fasta_file']), True),
                ('results/variants_all_clean.fa.gz', True)]
        if cfg.get('var_cdna_file'):
            out.append((str(cfg['var_cdna_file']), True))
    if cfg.get('add_existing_haplo'):
        out.append(('results/haplo_added_clean.fa.gz', True))
    out.append(('results/ref_contam_vcf_haplo_all_clean.fa.gz', True))
    out.append(('results/ref_contam_vcf_haplo_all_nodupl.fa.gz', True))
    final = str(cfg['final_fasta_file'])
    out.append((final, True))
    if cfg.get('fasta_simplify_headers'):
        if final.endswith('.gz'):
            out += [(strip_ext(final, 2) + '_simplified.fasta.gz', True), (strip_ext(final, 2) + '_header.tsv.gz', True)]
        else:
            out += [(strip_ext(final, 1) + '_simplified.fasta', True), (strip_ext(final, 1) + '_header.tsv', True)]
    return out


def write_report(path, results, times):
    with open(path, 'w') as fh:
        for arm in ARMS:
            if arm in times:
                fh.write('# %s: %.0f s wall, %.0f CPU-s, exit %d\n' % ((arm,) + times[arm]))
        fh.write('file\tlegacy_sha256\tfast_sha256\tresult\tfirst_difference\n')
        for r in results:
            fh.write('\t'.join([r.name, r.sha_a or '-', r.sha_b or '-', r.status, r.difference() or '-']) + '\n')


def print_report(results, times):
    width = max(len(r.name) for r in results)
    print()
    row = '%-*s  %-12s  %-12s  %-9s  %s'
    print(row % (width, 'file', 'legacy', 'fast', 'result', 'first difference'))
    for r in results:
        print((row % (width, r.name, (r.sha_a or '-')[:12], (r.sha_b or '-')[:12], r.status, r.difference())).rstrip())
    for r in results:
        if r.line is not None:
            print('\n%s, line %d:\n  legacy: %s\n  fast:   %s'
                  % (r.name, r.line, '<end of file>' if r.text_a is None else r.text_a,
                     '<end of file>' if r.text_b is None else r.text_b))
    print()
    for arm in ARMS:
        if arm in times:
            wall, cpu, code = times[arm]
            say('%-6s %8.0f s wall %10.0f CPU-s%s' % (arm, wall, cpu, '' if code == 0 else '   FAILED (exit %d)' % code))


# ------------------------------------ running ------------------------------------

def child_env():
    """The environment of every command verify starts: this Python's bin directory first, so that snakemake and
    python3 come from the same environment; a fixed hash seed and no user site-packages."""
    env = dict(os.environ)
    bindir = os.path.dirname(sys.executable)
    env['PATH'] = bindir + os.pathsep + env.get('PATH', '')
    env['PYTHONHASHSEED'] = '0'
    env['PYTHONNOUSERSITE'] = '1'
    return env


def snakemake(workdir, cores, log_path, extra=(), snakefile=SNAKEFILE):
    # the Snakefile reads config.yaml from the work directory
    cmd = [sys.executable, '-m', 'snakemake', '--snakefile', snakefile, '--directory', workdir,
           '--cores', str(cores)] + list(extra)
    with open(log_path, 'w') as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=child_env())
    return proc.returncode


def timed_snakemake(*args, **kwargs):
    """snakemake() plus its wall time and the CPU time of everything it ran, in seconds."""
    before, start = resource.getrusage(resource.RUSAGE_CHILDREN), time.monotonic()
    code = snakemake(*args, **kwargs)
    after, wall = resource.getrusage(resource.RUSAGE_CHILDREN), time.monotonic() - start
    cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
    return code, wall, cpu


def fast_src_workflow(run, fast_src):
    """A copy of the workflow whose src/ is a copy of fast_src (for --fast-src): legacy/, tools/ and envs/ link to
    this checkout's. src/ is copied, not linked: the fast scripts find legacy/src and tools/ next to their own
    directory, after resolving links."""
    wf = os.path.join(run.dir, 'workflow_fast')
    os.makedirs(wf)
    shutil.copy(SNAKEFILE, os.path.join(wf, 'Snakefile'))
    for name in ('legacy', 'tools', 'envs'):
        os.symlink(os.path.join(ROOT, name), os.path.join(wf, name))
    shutil.copytree(fast_src, os.path.join(wf, 'src'), ignore=shutil.ignore_patterns('__pycache__'))
    return os.path.join(wf, 'Snakefile')


def tail(path, n=30):
    with open(path, errors='replace') as fh:
        return ''.join(fh.readlines()[-n:])


FALLBACK = 'ProHap-fast: falling back to the reference code'     # src/modules/fallback.py


def fallbacks(path):
    """The lines of a run's log in which a fast script handed its input to the reference script."""
    with open(path, errors='replace') as fh:
        return [line.rstrip('\n') for line in fh if FALLBACK in line]


def prepare_list_with_rule(run):
    """Build a transcript list that a rule makes in the work directory, once, with the legacy arm's config."""
    def prepare(target):
        prep = os.path.join(run.dir, 'prep')
        os.makedirs(prep, exist_ok=True)
        write_yaml(os.path.join(prep, 'config.yaml'), run.arm_config('legacy'))
        log = os.path.join(run.dir, 'prep.log')
        say('building the transcript list %s once, for --transcripts' % target)
        if snakemake(prep, run.args.cores, log, ['--notemp', target]) != 0:
            raise VerifyError('building %s failed; see %s:\n%s' % (target, log, tail(log)))
        return os.path.join(prep, target)
    return prepare


def main(argv=None):
    args = parse_args(argv)
    cwd = os.getcwd()
    try:
        cfg = load_config(os.path.join(cwd, 'config.yaml'))
        run = Run(cfg, args, cwd)
        os.makedirs(run.input_dir)
        say('run directory', os.path.relpath(run.dir, cwd))
        say('chromosomes', ', '.join(str(c) for c in run.chromosomes))

        try:
            names = set(read_metadata(str(cfg['sample_metadata_file']))['Sample name'])
        except Exception:  # the workflow reports an unreadable metadata file itself
            names = None
        warnings, notes = preflight.check(run.chunks, names, shutil.which('tabix', path=child_env()['PATH']))
        for w in warnings:
            say('warning:', w)
        for n in notes:
            say('note:', n)
        if cfg.get('restrict_to_exons'):
            say('note: restrict_to_exons is on, so both arms read only the records that overlap an exon of a listed '
                'transcript. Records whose ID is "." are named by their line number in the filtered VCF, so their '
                'IDs differ from those of a run without the option (not between the arms). Set IDs first '
                '(bcftools annotate --set-id) if results must match a run without it.')

        samples = select_samples(run)
        if samples is not None:
            subset_vcfs(run, samples)
            rows = subset_metadata(run, samples)
            say('samples: %d (%s), %d of them in sample_metadata_file'
                % (len(samples), 'seed %d' % args.seed if args.samples else args.sample_list, rows))
        if args.transcripts:
            counts = subset_transcripts(run, args.transcripts, prepare_list_with_rule(run))
            say('transcripts per chromosome: ' + ', '.join('%s: %d' % kv for kv in counts.items()))

        run.write_configs()
        snakefiles = {'legacy': SNAKEFILE, 'fast': fast_src_workflow(run, args.fast_src) if args.fast_src else SNAKEFILE}
        if args.fast_src:
            say('the fast arm runs the scripts in', args.fast_src)

        if args.dry_run:
            failed = False
            for arm in ARMS:
                log = os.path.join(run.dir, arm + '.dryrun.log')
                if snakemake(run.arm_dir(arm), args.cores, log, ['--dry-run', '--quiet'], snakefiles[arm]) != 0:
                    say('%s arm: dry run FAILED; see %s:\n%s' % (arm, os.path.relpath(log, cwd), tail(log)))
                    failed = True
                else:
                    say('%s arm: dry run OK' % arm)
            return 2 if failed else 0

        if samples is None and not args.transcripts:
            first = next((paths[0] for paths in run.chunks.values() if paths), None)
            n = len((preflight.read_samples(first) if first else None) or [])
            say('warning: the legacy arm runs at full scale: %d samples, every transcript of %d chromosome(s). The '
                'reference takes about 2,000 CPU-s per 60 transcripts at 26k samples, and more than linearly longer '
                'with more samples. For a first check use --samples 2000 --transcripts 50, or --transcripts 20 to '
                'keep every sample.' % (n, len(run.chromosomes)))

        times = {}
        for arm in ARMS:
            log = os.path.join(run.dir, arm + '.log')
            say('running the %s arm; log: %s' % (arm, os.path.relpath(log, cwd)))
            code, wall, cpu = timed_snakemake(run.arm_dir(arm), args.cores, log, ['--notemp'], snakefiles[arm])
            times[arm] = (wall, cpu, code)
            handed_over = fallbacks(log) if arm == 'fast' else []
            if handed_over:
                say('note: %d step(s) of the fast arm handed their input to the reference code, which then wrote '
                    'their outputs (input the fast code does not handle, IMPLEMENTATION_NOTES section 3.7):\n  %s'
                    % (len(handed_over), '\n  '.join(line[:300] for line in handed_over)))
            if code != 0:
                say('%s arm FAILED (exit %d); the end of its log:\n%s' % (arm, code, tail(log)))
                if arm == 'legacy':
                    say('the reference failed, so there is nothing to compare against; kept', os.path.relpath(run.dir, cwd))
                    return 2

        results = [compare.compare(name, os.path.join(run.arm_dir('legacy'), name),
                                   os.path.join(run.arm_dir('fast'), name), required)
                   for name, required in expected_outputs(run.arm_config('fast'), run.refs)]
        results = [r for r in results if r.status != compare.ABSENT]
        report = os.path.join(run.dir, 'report.tsv')
        write_report(report, results, times)
        print_report(results, times)
        identical = all(r.ok for r in results) and times['fast'][2] == 0
        different = [r for r in results if not r.ok]
        if identical:
            say('IDENTICAL: all %d outputs' % len(results))
            if not args.keep:
                for name in os.listdir(run.dir):
                    if name != 'report.tsv':
                        path = os.path.join(run.dir, name)
                        shutil.rmtree(path) if os.path.isdir(path) and not os.path.islink(path) else os.remove(path)
                say('removed the work directories (--keep keeps them); report:', os.path.relpath(report, cwd))
            else:
                say('kept', os.path.relpath(run.dir, cwd))
            return 0
        if times['fast'][2] != 0:
            say('FAILED: the fast arm failed, and %d of %d outputs differ or are missing; kept %s'
                % (len(different), len(results), os.path.relpath(run.dir, cwd)))
        else:
            say('DIFFERENT: %d of %d outputs; kept %s' % (len(different), len(results), os.path.relpath(run.dir, cwd)))
        return 1
    except VerifyError as e:
        say('error:', e)
        return 2


if __name__ == '__main__':
    sys.exit(main())
