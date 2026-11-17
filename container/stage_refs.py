#!/usr/bin/env python3
"""Stage bundled Ensembl reference files into the ProHap work directory.

The Snakefile obtains its reference files through `download_*` rules that call
wget. In an air-gapped environment those rules must never run, so this script
places their outputs on disk beforehand; Snakemake then treats them as already
satisfied.

Staging is lazy. Each reference is only decompressed if the derived artifact
that consumes it is missing from the work directory, so a warm cache costs
nothing. This matters because the uncompressed references total roughly 2 GB.

Reads the ProHap config to learn the annotation filename, the Ensembl release
and the chromosome list, so that the staged paths match what the rules expect.
"""

import argparse
import gzip
import os
import shutil
import sys

import yaml


def log(msg):
    print('[stage_refs] ' + msg, flush=True)


def fail(msg):
    print('[stage_refs] ERROR: ' + msg, file=sys.stderr, flush=True)
    sys.exit(1)


def inflate(src_gz, dest):
    """Decompress src_gz to dest, atomically."""
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + '.partial'
    with gzip.open(src_gz, 'rb') as fin, open(tmp, 'wb') as fout:
        shutil.copyfileobj(fin, fout, length=8 * 1024 * 1024)
    os.replace(tmp, dest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-c', '--config', required=True, help='ProHap config.yaml')
    parser.add_argument('-r', '--refs', required=True, help='directory of gzipped Ensembl reference files')
    parser.add_argument('-w', '--workdir', required=True, help='ProHap working directory')
    args = parser.parse_args()

    with open(args.config) as fh:
        config = yaml.safe_load(fh)

    for key in ('annotationFilename', 'ensembl_release', 'chromosomes'):
        if key not in config:
            fail('config is missing the required key "%s"' % key)

    annot = config['annotationFilename']
    release = str(config['ensembl_release'])
    chromosomes = [str(c) for c in config['chromosomes']]

    # The single most common way to break an air-gapped run: the bundled
    # references are one release and the config asks for another, so the
    # download rules fire and fail with no network.
    if release not in annot:
        log('WARNING: ensembl_release is %s but annotationFilename is "%s" -'
            ' these usually have to agree.' % (release, annot))

    work = args.workdir
    gtf = os.path.join(work, 'data', 'gtf', annot + '.gtf')

    # Map each bundled reference to the derived artifacts that consume it.
    # If every consumer already exists, the reference is not needed at all.
    whole_db = os.path.join(work, 'data', 'gtf', annot + '.db')
    chr_dbs = [os.path.join(work, 'data', 'gtf', '%s_chr%s.db' % (annot, c)) for c in chromosomes]

    plan = [
        # (bundled file, staged destination, [consumers])
        (annot + '.gtf.gz', gtf, [whole_db] + chr_dbs),
        ('Homo_sapiens.GRCh38.cdna.all.fa.gz',
         os.path.join(work, 'data', 'fasta', 'Homo_sapiens.GRCh38.cdna.all.fa'),
         [os.path.join(work, 'data', 'fasta', 'total_cdnas_%s.fa' % release)]),
        ('Homo_sapiens.GRCh38.ncrna.fa.gz',
         os.path.join(work, 'data', 'fasta', 'Homo_sapiens.GRCh38.ncrna.fa'),
         [os.path.join(work, 'data', 'fasta', 'total_cdnas_%s.fa' % release)]),
        ('Homo_sapiens.GRCh38.pep.all.fa.gz',
         os.path.join(work, 'data', 'fasta', 'Homo_sapiens.GRCh38.pep.all.fa'),
         [os.path.join(work, 'data', 'fasta',
                       'ensembl_reference_proteinDB_%s_tagged.fa' % release)]),
    ]

    staged = 0
    for bundled, dest, consumers in plan:
        if os.path.exists(dest):
            continue
        missing = [c for c in consumers if not os.path.exists(c)]
        if not missing:
            log('skipping %s - all derived files already cached' % os.path.basename(dest))
            continue

        src = os.path.join(args.refs, bundled)
        if not os.path.exists(src):
            fail('need %s but %s is not in the image.\n'
                 '  Missing derived file(s): %s\n'
                 '  Either rebuild the image with a matching ENSEMBL_RELEASE, or\n'
                 '  point PROHAP_REFS at a directory containing the gzipped file.'
                 % (os.path.basename(dest), src, ', '.join(missing[:3])))

        log('staging %s -> %s' % (bundled, dest))
        inflate(src, dest)
        staged += 1

    if staged == 0:
        log('cache is warm, no references needed staging')


if __name__ == '__main__':
    main()
