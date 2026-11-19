#!/usr/bin/env python3
"""Time the filter and prohap.py of each implementation on simulated data: wall time, CPU time, memory, and the
phases of prohap.py.

    python3 bench/bench.py --samples 2500,5000 --arms legacy,fast
    python3 bench/bench.py --data DIR --arms legacy --threads 4

For every dataset, chromosome and arm, the harness runs the two steps of the workflow that scale with the cohort,
with the arguments the Snakefile gives them: vcf_filter_fix.py (the `filter` step), then prohap.py (`prohap`). Each
runs with PYTHONHASHSEED=0 (through tools/run_fork.py on macOS), in its own work directory. Per step it records:

    wall      seconds from start to exit
    cpu       user + system seconds of the step and every process it waited for (its Pool workers), from the
              rusage that wait4() returns for the step's process
    tree_mb   peak of the summed resident memory of the step's process tree, sampled every 0.25 s with ps (forked
              workers share pages, so this overstates; it is what a job scheduler sees)
    max_mb    the largest resident memory of any single process (ru_maxrss)
    phases    for prohap.py, wall seconds between the lines it prints: setup, sweep (assigning variants to
              transcripts), haplotypes (co-occurrence and frequencies), cdna (reading the cDNA FASTA), proteins,
              write. The fast implementation reads the cDNA during setup and computes haplotypes and proteins while
              it reads the VCF, so its phases are setup, cdna, sweep (everything streamed) and write.

When two arms ran, their outputs (the filtered VCF, haplotype table, FASTA, cDNA FASTA and removed-samples log) are
compared after decompression; the identity check proper is tests/legacy_vs_fast.sh.

Data: --data DIR (any directory written by bench/simulate.py; every chromosome in its transcripts.csv), or
--samples N,... for the benchmark block: --transcripts consecutive MANE transcripts of --chrom from --first (default
60 on chromosome 22 from ENST00000402142, 217 kb of exons, about 1/315 of the autosomal exome), simulated once into
bench/data/ and reused while bench/simulate.py is unchanged.

Results are printed as a table and appended to --tsv (default bench/results.tsv).
"""
import argparse
import datetime
import gzip
import hashlib
import json
import os
import platform
import resource
import shutil
import subprocess
import sys
import tempfile
import threading
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SIMULATE = os.path.join(ROOT, 'bench', 'simulate.py')
ARMS = {'legacy': os.path.join(ROOT, 'legacy', 'src'), 'fast': os.path.join(ROOT, 'src')}
PHASE_MARKERS = [               # (phase that starts, text of the line prohap.py prints when it starts)
    ('sweep', 'Assigning variants to transcripts.'),
    ('haplotypes', 'Computing the co-occurrence of alleles.'),
    ('cdna', 'Reading {cdna}'),
    ('proteins', 'Creating haplotype database.'),
    ('write', 'Storing the result metadata'),
]
COLUMNS = ['date', 'host', 'dataset', 'chrom', 'samples', 'transcripts', 'records', 'arm', 'step', 'threads', 'wall',
           'cpu', 'tree_mb', 'max_mb', 'phases']
SAMPLE_EVERY = 0.25


def default_refs():
    root = os.environ.get('PROHAP_REFS', '')
    found = sorted(d for d in (os.listdir(root) if root and os.path.isdir(root) else []) if d.startswith('ensembl.v'))
    return os.path.join(root, found[0]) if len(found) == 1 else None


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def sha256_content(path):
    """sha256 of the decompressed bytes (gzip by magic number)."""
    with open(path, 'rb') as fh:
        gz = fh.read(2) == b'\x1f\x8b'
    h = hashlib.sha256()
    with (gzip.open(path, 'rb') if gz else open(path, 'rb')) as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------- data

def standard_dataset(args, samples):
    """The benchmark block at this sample count, simulated into bench/data/ unless already there."""
    name = 'chr%s-%s-%d-n%d-%s-seed%d' % (args.chrom, args.first, args.transcripts, samples, args.ids, args.seed)
    path = os.path.join(ROOT, 'bench', 'data', name)
    summary = os.path.join(path, 'chr%s.json' % args.chrom)
    try:
        with open(summary) as fh:
            if json.load(fh)['simulate_sha256'] == sha256_file(SIMULATE):
                return path
    except (OSError, ValueError, KeyError):
        pass
    shutil.rmtree(path, ignore_errors=True)
    print('simulating %s ...' % name, flush=True)
    subprocess.run([sys.executable, SIMULATE, '--refs', args.refs, '--chrom', args.chrom, '--first', args.first,
                    '--transcripts', str(args.transcripts), '--samples', str(samples), '--ids', args.ids,
                    '--seed', str(args.seed), '--out', path], check=True)
    return path


def dataset_chromosomes(path):
    with open(os.path.join(path, 'transcripts.csv')) as fh:
        chroms = [line.split(',', 1)[0] for line in list(fh)[1:]]
    return list(dict.fromkeys(chroms))


def dataset_info(path, chrom):
    try:
        with open(os.path.join(path, 'chr%s.json' % chrom)) as fh:
            info = json.load(fh)
        samples = int(info['command'].split('--samples ')[1].split()[0])
        return samples, len(info['block']), info['records']
    except (OSError, ValueError, KeyError, IndexError):
        return '', '', ''


# ---------------------------------------------------------------- measuring one step

def tree_rss_kb(root_pid):
    """Summed RSS (kB) of root_pid and all its descendants."""
    try:
        out = subprocess.run(['ps', '-A', '-o', 'pid=,ppid=,rss='], capture_output=True, text=True).stdout
    except OSError:
        return 0
    children = {}
    rss = {}
    for line in out.splitlines():
        parts = line.split()
        if len(parts) == 3:
            pid, ppid, kb = map(int, parts)
            children.setdefault(ppid, []).append(pid)
            rss[pid] = kb
    total = 0
    stack = [root_pid]
    while stack:
        pid = stack.pop()
        total += rss.get(pid, 0)
        stack.extend(children.get(pid, []))
    return total


def run_step(cmd, cwd, env, log_path):
    """Run cmd; returns (exit code, wall s, cpu s, peak tree RSS MB, max RSS MB, [(seconds, stdout line)])."""
    start = time.monotonic()
    proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                            errors='replace')
    peak = [0]
    done = threading.Event()

    def sample():
        while not done.is_set():
            peak[0] = max(peak[0], tree_rss_kb(proc.pid))
            done.wait(SAMPLE_EVERY)

    sampler = threading.Thread(target=sample, daemon=True)
    sampler.start()
    lines = []
    with open(log_path, 'w') as log:
        for line in proc.stdout:
            lines.append((time.monotonic() - start, line.rstrip('\n')))
            log.write(line)
    # wait4 returns the rusage of this child and the processes it waited for (prohap.py's Pool workers), and not
    # that of the ps calls above, which RUSAGE_CHILDREN of this process would include
    _, status, usage = os.wait4(proc.pid, 0)
    proc.returncode = os.waitstatus_to_exitcode(status)
    wall = time.monotonic() - start
    done.set()
    sampler.join()
    max_rss = usage.ru_maxrss / (1 << 20 if sys.platform == 'darwin' else 1 << 10)     # bytes on macOS, kB on Linux
    return proc.returncode, wall, usage.ru_utime + usage.ru_stime, peak[0] / 1024, max_rss, lines


def phases(lines, wall, cdna):
    marks = [(0.0, 'setup')]
    for name, text in PHASE_MARKERS:
        text = text.format(cdna=cdna)
        hit = next((t for t, line in lines if text in line), None)
        if hit is not None:
            marks.append((hit, name))
    marks.sort()
    marks.append((wall, None))
    return ';'.join('%s=%.1f' % (marks[i][1], marks[i + 1][0] - marks[i][0]) for i in range(len(marks) - 1))


# ---------------------------------------------------------------- main

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], epilog=__doc__.split('\n\n', 1)[1],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--refs', default=default_refs(),
                        help='references directory ensembl.v<release> (default: the one in $PROHAP_REFS)')
    data = parser.add_mutually_exclusive_group(required=True)
    data.add_argument('--data', action='append', help='a dataset directory written by bench/simulate.py (repeatable)')
    data.add_argument('--samples', help='comma-separated sample counts for the benchmark block')
    parser.add_argument('--chrom', default='22', help='benchmark block: chromosome (default 22)')
    parser.add_argument('--first', default='ENST00000402142', help='benchmark block: first transcript')
    parser.add_argument('--transcripts', type=int, default=60, help='benchmark block: transcripts (default 60)')
    parser.add_argument('--ids', default='dot', help='benchmark block: simulate.py --ids (default dot)')
    parser.add_argument('--seed', type=int, default=0, help='benchmark block: simulate.py --seed (default 0)')
    parser.add_argument('--arms', default='legacy,fast',
                        help='comma-separated: legacy (legacy/src), fast (src/) or a script directory '
                             '(default legacy,fast)')
    parser.add_argument('--threads', type=int, default=1, help='prohap.py -threads (default 1)')
    parser.add_argument('--af', default='0', help='vcf_filter_fix.py -af (default 0: keep every variant)')
    parser.add_argument('--min-freq', default='0', help='prohap.py -min_hap_freq (default 0)')
    parser.add_argument('--min-count', default='0', help='prohap.py -min_hap_count (default 0)')
    parser.add_argument('--cdna', action='store_true', help='also write the cDNA FASTA')
    parser.add_argument('--work', help='work directory (default: a new temporary directory, deleted at the end '
                                       'unless --keep)')
    parser.add_argument('--keep', action='store_true', help='keep the work directory')
    parser.add_argument('--tsv', default=os.path.join(ROOT, 'bench', 'results.tsv'),
                        help='append results here (default bench/results.tsv)')
    args = parser.parse_args(argv)
    if not args.refs:
        parser.error('--refs is required ($PROHAP_REFS does not hold exactly one ensembl.v<release>)')
    args.refs = os.path.abspath(os.path.expanduser(args.refs))
    with open(os.path.join(args.refs, 'release.txt')) as fh:
        release = fh.read().strip()
    cdna = os.path.join(args.refs, 'total_cdnas_%s.fa' % release)

    arms = []
    for arm in args.arms.split(','):
        src = ARMS.get(arm, os.path.abspath(arm))
        if not os.path.isfile(os.path.join(src, 'prohap.py')):
            parser.error('arm %s: no prohap.py in %s' % (arm, src))
        arms.append((arm if arm in ARMS else os.path.basename(src.rstrip('/')), src))
    datasets = [os.path.abspath(d) for d in args.data] if args.data else \
        [standard_dataset(args, int(n)) for n in args.samples.split(',')]

    work = os.path.abspath(args.work) if args.work else tempfile.mkdtemp(prefix='prohap-bench.')
    os.makedirs(work, exist_ok=True)
    env = dict(os.environ, PYTHONHASHSEED='0', PYTHONNOUSERSITE='1', PYTHONUNBUFFERED='1')
    runner = [sys.executable] + ([os.path.join(ROOT, 'tools', 'run_fork.py')] if sys.platform == 'darwin' else [])
    host = platform.node().split('.')[0] + ' ' + platform.machine()
    date = datetime.datetime.now().isoformat(timespec='seconds')
    rows = []
    failed = False
    comparisons = []

    print('work directory:', work)
    for path in datasets:
        dataset = os.path.basename(path.rstrip('/'))
        for chrom in dataset_chromosomes(path):
            samples, n_transcripts, records = dataset_info(path, chrom)
            vcfs = sorted(f for f in os.listdir(path) if f.endswith('.vcf.gz') and
                          (f == 'chr%s.vcf.gz' % chrom or f.startswith('chr%s.' % chrom)))
            if not vcfs:
                print('%s: no VCF for chromosome %s, skipped' % (dataset, chrom))
                continue
            outputs = {}
            for arm, src in arms:
                d = os.path.join(work, dataset, 'chr' + chrom, arm)
                shutil.rmtree(d, ignore_errors=True)
                os.makedirs(os.path.join(d, 'tmp'))
                glob = os.path.join(path, 'chr%s.vcf.gz' % chrom) if vcfs == ['chr%s.vcf.gz' % chrom] \
                    else os.path.join(path, 'chr%s.*vcf.gz' % chrom)
                steps = [
                    ('filter', runner + [os.path.join(src, 'vcf_filter_fix.py'), '-i', glob, '-chr', chrom,
                                         '-af', args.af, '-af_field', 'AF', '-o', 'filtered.vcf.gz']),
                    ('prohap', runner + [os.path.join(src, 'prohap.py'), '-i', 'filtered.vcf.gz', '-db',
                                         os.path.join(args.refs, 'Homo_sapiens.GRCh38.%s.chr_patch_hapl_scaff_chr%s.db'
                                                      % (release, chrom)),
                                         '-transcripts', os.path.join(path, 'transcripts.csv'), '-cdna', cdna,
                                         '-s', os.path.join(path, 'samples.tsv'), '-chr', chrom,
                                         '-min_hap_freq', args.min_freq, '-min_hap_count', args.min_count,
                                         '-acc_prefix', 'enshap_' + chrom, '-id_prefix', 'haplo_chr' + chrom,
                                         '-require_start', '1', '-ignore_UTR', '1', '-skip_start_lost', '1',
                                         '-x_par1_to', '2781479', '-x_par2_from', '155701383',
                                         '-threads', str(args.threads), '-log', 'prohap_chr%s.log' % chrom,
                                         '-tmp_dir', 'tmp', '-output_csv', 'haplo_chr%s.tsv.gz' % chrom,
                                         '-output_fasta', 'haplo_chr%s.fa' % chrom]
                               + (['-output_cdna_fasta', 'haplo_cdna_chr%s.fa' % chrom] if args.cdna else [])),
                ]
                for step, cmd in steps:
                    print('%s chr%s %s %s ...' % (dataset, chrom, arm, step), end=' ', flush=True)
                    code, wall, cpu, tree, maxrss, lines = run_step(cmd, d, env, os.path.join(d, step + '.out'))
                    print('%.1f s wall, %.1f s CPU, %.0f MB%s' % (wall, cpu, tree,
                                                                  '' if code == 0 else ', FAILED (exit %d)' % code),
                          flush=True)
                    rows.append([date, host, dataset, chrom, samples, n_transcripts, records, arm, step, args.threads,
                                 '%.1f' % wall, '%.1f' % cpu, '%.0f' % tree, '%.0f' % maxrss,
                                 phases(lines, wall, cdna) if step == 'prohap' else ''])
                    if code != 0:
                        failed = True
                        print('    see', os.path.join(d, step + '.out'))
                        break
                else:
                    outputs[arm] = {f: sha256_content(os.path.join(d, f)) for f in sorted(os.listdir(d))
                                    if f not in ('tmp', 'filter.out', 'prohap.out')}
            if len(outputs) > 1:
                names = list(outputs)
                files = sorted(set().union(*outputs.values()))
                same = all(len({o.get(f) for o in outputs.values()}) == 1 for f in files)
                comparisons.append((dataset, chrom, ' vs '.join(names), same,
                                    [f for f in files if len({o.get(f) for o in outputs.values()}) > 1]))

    widths = [max(len(str(r[i])) for r in rows + [COLUMNS]) for i in range(len(COLUMNS))]
    print()
    for r in [COLUMNS] + rows:
        print('  '.join(str(v).ljust(widths[i]) for i, v in enumerate(r) if COLUMNS[i] not in ('date', 'host')))
    for dataset, chrom, names, same, diff in comparisons:
        print('%s chr%s, %s: %s' % (dataset, chrom, names, 'outputs IDENTICAL' if same else
                                    'outputs DIFFERENT: ' + ', '.join(diff)))
    if args.tsv:
        new = not os.path.exists(args.tsv)
        with open(args.tsv, 'a') as fh:
            if new:
                fh.write('\t'.join(COLUMNS) + '\n')
            for r in rows:
                fh.write('\t'.join(map(str, r)) + '\n')
        print('appended to', args.tsv)
    if args.keep or failed:
        print('kept', work)
    elif not args.work:
        shutil.rmtree(work)
    return 1 if failed or any(not c[3] for c in comparisons) else 0


if __name__ == '__main__':
    sys.exit(main())
