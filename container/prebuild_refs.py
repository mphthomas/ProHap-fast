#!/usr/bin/env python3
"""Build every reference ProHap derives from an Ensembl release, once, into a read-only directory.

Input: the four Ensembl downloads of one release (the GTF `chr_patch_hapl_scaff`, cDNA, ncRNA and protein FASTA),
gzipped, as the Ensembl FTP serves them, plus the contaminants FASTA. Output: `<out>/ensembl.v<release>/` with

    Homo_sapiens.GRCh38.<rel>.chr_patch_hapl_scaff.db           whole-genome gffutils DB
    Homo_sapiens.GRCh38.<rel>.chr_patch_hapl_scaff_chr<c>.db    one per chromosome
    transcripts_reference_<rel>.csv                             MANE Select list (Ensembl canonical without MANE)
    total_cdnas_<rel>.fa                                        ncRNA then cDNA
    ensembl_reference_proteinDB_<rel>_tagged.fa                 reference proteome, tagged
    crap_tagged.fa.gz                                           contaminants, tagged
    release.txt                                                 the release number
    SOURCES.sha256                                              sha256 of the inputs these were built from
    MANIFEST.sha256                                             sha256 of every file above (`sha256sum -c`)

Every file is built by the frozen reference scripts in legacy/src, called with the arguments the upstream rules
used (`parse_gtf_whole`, `split_gtf` + `parse_gtf_chromosome`, `merge_cdnas_fasta`, `reference_filter_format`,
`default_transcript_list` with only_MANE_select, `contaminants_fix_headers`), so the workflow reads exactly what it
used to build in each work directory. The directory appears only once it is complete: it is built under
`ensembl.v<rel>.partial/` and renamed at the end.

    python3 container/prebuild_refs.py --release 113 --downloads DIR --out REFS_ROOT [-j 8]
"""
import argparse
import concurrent.futures
import gzip
import hashlib
import os
import re
import shlex
import shutil
import sqlite3
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.realpath(__file__))
REPO = os.path.dirname(HERE)
DEFAULT_SRC = os.path.join(REPO, 'legacy', 'src')
DEFAULT_CONTAMINANTS = os.path.join(REPO, 'crap.fasta')
# Every chromosome a config can name: GTF seqnames of the primary assembly.
DEFAULT_CHROMOSOMES = [str(c) for c in range(1, 23)] + ['X', 'Y', 'MT']


def annotation_name(release):
    return 'Homo_sapiens.GRCh38.%s.chr_patch_hapl_scaff' % release


def download_names(release):
    """The Ensembl FTP file names, keyed by what each one is."""
    return {
        'gtf': annotation_name(release) + '.gtf.gz',
        'cdna': 'Homo_sapiens.GRCh38.cdna.all.fa.gz',
        'ncrna': 'Homo_sapiens.GRCh38.ncrna.fa.gz',
        'pep': 'Homo_sapiens.GRCh38.pep.all.fa.gz',
    }


def output_names(release, chromosomes):
    """File names inside ensembl.v<release>/, the contract with the Snakefile."""
    annot = annotation_name(release)
    names = {
        'whole_db': annot + '.db',
        'mane': 'transcripts_reference_%s.csv' % release,
        'cdna': 'total_cdnas_%s.fa' % release,
        'proteome': 'ensembl_reference_proteinDB_%s_tagged.fa' % release,
        'contaminants': 'crap_tagged.fa.gz',
    }
    for c in chromosomes:
        names['chr_db:' + c] = '%s_chr%s.db' % (annot, c)
    return names


def log(msg):
    print('[prebuild_refs] ' + msg, flush=True)


def fail(msg):
    print('[prebuild_refs] ERROR: ' + msg, file=sys.stderr, flush=True)
    sys.exit(1)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def gunzip(src, dest):
    """What `gunzip` did to each download in the upstream rules."""
    with gzip.open(src, 'rb') as fin, open(dest, 'wb') as fout:
        shutil.copyfileobj(fin, fout, length=8 * 1024 * 1024)


class Builder:
    def __init__(self, args):
        self.args = args
        self.python = args.python
        self.src = args.src
        self.work = args.work
        self.logs = os.path.join(args.work, 'logs')
        os.makedirs(self.logs, exist_ok=True)
        # The scripts run in the environment the workflow's rules see (Snakefile shell.prefix).
        self.env = dict(os.environ, PYTHONHASHSEED='0', PYTHONNOUSERSITE='1')

    def run(self, name, cmd):
        """Run one shell command (bash strict mode, like a Snakemake rule); its output goes to logs/<name>.log."""
        start = time.monotonic()
        log_path = os.path.join(self.logs, name + '.log')
        with open(log_path, 'w') as log_fh:
            log_fh.write('$ ' + cmd + '\n')
            log_fh.flush()
            proc = subprocess.run(['bash', '-c', 'set -euo pipefail; ' + cmd], stdout=log_fh,
                                  stderr=subprocess.STDOUT, env=self.env, cwd=self.work)
        took = time.monotonic() - start
        if proc.returncode != 0:
            with open(log_path) as fh:
                tail = fh.readlines()[-20:]
            raise RuntimeError('%s failed (exit %d) after %.1f s; last lines of %s:\n%s'
                               % (name, proc.returncode, took, log_path, ''.join(tail)))
        log('%-28s %7.1f s' % (name, took))
        return took

    def script(self, name, *args):
        return ' '.join([shlex.quote(self.python), shlex.quote(os.path.join(self.src, name))]
                        + [shlex.quote(a) for a in args])


def check_db(path):
    """The workflow opens the DBs from a read-only directory: they must open read-only and not be in WAL mode."""
    con = sqlite3.connect('file:' + path + '?mode=ro', uri=True)
    try:
        mode = con.execute('PRAGMA journal_mode').fetchone()[0]
        n = con.execute('SELECT count(*) FROM features').fetchone()[0]
    finally:
        con.close()
    if mode.lower() == 'wal':
        raise RuntimeError('%s is in WAL journal mode, which needs write access to read' % path)
    if n == 0:
        raise RuntimeError('%s has no features' % path)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--release', required=True, help='Ensembl release, e.g. 113')
    parser.add_argument('--downloads', required=True, help='directory holding the four gzipped Ensembl downloads')
    parser.add_argument('--out', required=True, help='refs root; the result is <out>/ensembl.v<release>/')
    parser.add_argument('--contaminants', default=DEFAULT_CONTAMINANTS,
                        help='contaminants FASTA (default: crap.fasta in the repository)')
    parser.add_argument('--chromosomes', default=','.join(DEFAULT_CHROMOSOMES),
                        help='chromosomes that get their own DB, comma-separated (default: 1-22,X,Y,MT)')
    parser.add_argument('--src', default=DEFAULT_SRC, help='ProHap scripts to build with (default: legacy/src)')
    parser.add_argument('--python', default=sys.executable, help='interpreter for the ProHap scripts')
    parser.add_argument('--tmp', default=None,
                        help='scratch space for the decompressed downloads, about 3 GB (default: inside --out)')
    parser.add_argument('-j', '--jobs', type=int, default=min(8, os.cpu_count() or 1), help='parallel builds')
    parser.add_argument('--force', action='store_true', help='replace an existing ensembl.v<release>/')
    args = parser.parse_args()

    release = str(args.release)
    if not release.isdigit():
        fail('--release must be a number, not %r' % release)
    chromosomes = [c.strip() for c in args.chromosomes.split(',') if c.strip()]
    if len(set(chromosomes)) != len(chromosomes):
        fail('--chromosomes lists a chromosome twice')
    for c in chromosomes:
        if not re.fullmatch(r'[A-Za-z0-9_]+', c):
            fail('chromosome names are letters, digits and _ only, not %r' % c)
    args.src = os.path.realpath(args.src)
    args.python = os.path.abspath(shutil.which(args.python) or args.python)  # not realpath: keeps a venv

    downloads = {k: os.path.join(args.downloads, v) for k, v in download_names(release).items()}
    for path in list(downloads.values()) + [args.contaminants]:
        if not os.path.isfile(path):
            fail('missing input %s' % path)
    for name in ('parse_gtf.py', 'get_reference_ENST.py', 'fasta_format_headers.py'):
        if not os.path.isfile(os.path.join(args.src, name)):
            fail('%s is not in --src %s' % (name, args.src))

    os.makedirs(args.out, exist_ok=True)
    final = os.path.join(args.out, 'ensembl.v' + release)
    partial = final + '.partial'
    if os.path.exists(final):
        if not args.force:
            fail('%s already exists; pass --force to rebuild it' % final)
    if os.path.exists(partial):
        shutil.rmtree(partial)
    os.makedirs(partial)
    args.work = os.path.join(args.tmp or partial, '.prebuild_work' if args.tmp is None else 'prebuild_work.' + release)
    if os.path.exists(args.work):
        shutil.rmtree(args.work)
    os.makedirs(args.work)

    t0 = time.monotonic()
    log('Ensembl release %s, scripts %s, %d jobs' % (release, args.src, args.jobs))
    log('output %s' % final)

    # Input checksums first: they record exactly what the references were built from.
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        source_paths = dict(downloads, contaminants=args.contaminants)
        source_sums = dict(zip(source_paths, pool.map(sha256_file, source_paths.values())))

    # The upstream rules downloaded each file and ran gunzip on it; the scripts read the plain files.
    plain = {k: os.path.join(args.work, os.path.basename(v)[:-len('.gz')]) for k, v in downloads.items()}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        list(pool.map(lambda k: gunzip(downloads[k], plain[k]), plain))
    log('%-28s %7.1f s' % ('decompress downloads', time.monotonic() - t0))

    b = Builder(args)
    out = {k: os.path.join(partial, v) for k, v in output_names(release, chromosomes).items()}
    q = shlex.quote

    def chromosome_db(c):
        # rule split_gtf, then rule parse_gtf_chromosome
        chr_gtf = os.path.join(args.work, '%s_chr%s.gtf' % (annotation_name(release), c))
        took = b.run('split_gtf chr' + c,
                     'grep "^#" %s > %s; grep "^%s\\s" %s >> %s' % (q(plain['gtf']), q(chr_gtf), c,
                                                                     q(plain['gtf']), q(chr_gtf)))
        took += b.run('parse_gtf chr' + c, b.script('parse_gtf.py', '-i', chr_gtf, '-o', out['chr_db:' + c]))
        os.remove(chr_gtf)
        return took

    jobs = {
        # rule parse_gtf_whole: the longest job, so it starts first
        'whole_db': lambda: b.run('parse_gtf whole genome',
                                  b.script('parse_gtf.py', '-i', plain['gtf'], '-o', out['whole_db'])),
        # rule reference_filter_format
        'proteome': lambda: b.run('reference_filter_format',
                                  b.script('fasta_format_headers.py', '-i', plain['pep'], '-o', out['proteome'],
                                           '-t', '_ensref', '-use_ENST', '1')),
        # rule merge_cdnas_fasta
        'cdna': lambda: b.run('merge_cdnas_fasta',
                              'cat %s > %s; cat %s >> %s' % (q(plain['ncrna']), q(out['cdna']),
                                                             q(plain['cdna']), q(out['cdna']))),
        # rule contaminants_fix_headers
        'contaminants': lambda: b.run('contaminants_fix_headers',
                                      b.script('fasta_format_headers.py', '-i', args.contaminants,
                                               '-o', out['contaminants'], '-t', '_cont')),
    }
    for c in chromosomes:
        jobs['chr_db:' + c] = (lambda c=c: chromosome_db(c))

    step_time = 0.0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {name: pool.submit(fn) for name, fn in jobs.items()}
        try:
            # rule default_transcript_list with only_MANE_select: needs the whole DB and the tagged proteome
            futures['whole_db'].result()
            futures['proteome'].result()
            futures['mane'] = pool.submit(
                b.run, 'default_transcript_list',
                b.script('get_reference_ENST.py', '-i', out['proteome'], '-annot', out['whole_db'],
                         '-MANE', '1', '-o', out['mane']))
            for f in futures.values():
                step_time += f.result()
        except Exception as e:
            pool.shutdown(wait=False, cancel_futures=True)
            fail('%s\nThe partial build is kept in %s.' % (e, partial))

    shutil.rmtree(args.work)

    for k, path in out.items():
        if k == 'whole_db' or k.startswith('chr_db:'):
            check_db(path)

    with open(os.path.join(partial, 'release.txt'), 'w') as fh:
        fh.write(release + '\n')
    with open(os.path.join(partial, 'SOURCES.sha256'), 'w') as fh:
        for k in ('gtf', 'cdna', 'ncrna', 'pep', 'contaminants'):
            fh.write('%s  %s\n' % (source_sums[k], os.path.basename(source_paths[k])))

    names = sorted(os.listdir(partial))
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        sums = list(pool.map(lambda n: sha256_file(os.path.join(partial, n)), names))
    with open(os.path.join(partial, 'MANIFEST.sha256'), 'w') as fh:
        for name, digest in zip(names, sums):
            fh.write('%s  %s\n' % (digest, name))

    # Readable by any user: the image runs as whoever calls it.
    for name in os.listdir(partial):
        os.chmod(os.path.join(partial, name), 0o644)
    os.chmod(partial, 0o755)

    if os.path.exists(final):
        shutil.rmtree(final)
    os.rename(partial, final)

    size = sum(os.path.getsize(os.path.join(final, n)) for n in os.listdir(final))
    log('done in %.0f s (build steps add up to %.0f s), %d files, %.2f GB: %s'
        % (time.monotonic() - t0, step_time, len(os.listdir(final)), size / 1e9, final))


if __name__ == '__main__':
    main()
