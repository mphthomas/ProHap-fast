#!/usr/bin/env python3
"""Generate the data of a simulated fixture into tests/fixtures/sim/cache/<dataset>/, unless it is there already, and
print the directory.

    python3 tests/fixtures/sim/make_data.py DATASET [--refs DIR]
    python3 tests/fixtures/sim/make_data.py --list

The data is written by bench/simulate.py (VCFs, samples.tsv, transcripts.csv) or bench/simulate_fasta.py (dupseq)
from fixed seeds, so it is never committed: the cache is git-ignored and rebuilt whenever the dataset's definition
below, the simulator or the Ensembl release of the references changes (cache/<dataset>/stamp.json records them).
The fixtures and what each one tests are described in README.md next to this file.
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
CACHE = os.path.join(HERE, 'cache')
SIMULATE = os.path.join(ROOT, 'bench', 'simulate.py')
SIMULATE_FASTA = os.path.join(ROOT, 'bench', 'simulate_fasta.py')

# The benchmark block: 60 consecutive MANE transcripts of chromosome 22, 217 kb of exons (bench/bench.py).
BLOCK22 = ['--chrom', '22', '--first', 'ENST00000402142', '--transcripts', '60']

# dataset: samples in the VCFs, samples left out of the metadata (and as many metadata-only samples), and one
# bench/simulate.py argument list per chromosome; or a bench/simulate_fasta.py argument list ('fasta')
DATASETS = {
    # small tier: one chromosome per ID mode, about 30 kb of exons each; on X, 5 PAR1 transcripts (two of them
    # overlapping), one across the PAR1 boundary and 3 beyond it; chromosome 22 in two VCF chunks
    'mix': {'samples': 2500, 'unlisted': 5, 'chromosomes': [
        ['--chrom', '20', '--first', 'ENST00000381458', '--transcripts', '15', '--ids', 'ambiguous'],
        ['--chrom', '21', '--first', 'ENST00000397928', '--transcripts', '15', '--ids', 'rs'],
        ['--chrom', '22', '--first', 'ENST00000292778', '--transcripts', '15', '--ids', 'dot', '--chunks', '2'],
        ['--chrom', 'X', '--first', 'ENST00000313871', '--transcripts', '9', '--ids', 'dot', '--males', 'haploid'],
    ]},
    # medium tier: the benchmark block at growing cohort sizes (pandas reads the per-transcript tables in chunks of
    # 256, 128 and 64 rows), chrX with mirrored males, a merged table of more than 32,768 rows
    'n2500': {'samples': 2500, 'unlisted': 0, 'chromosomes': [BLOCK22 + ['--ids', 'dot']]},
    'n5000': {'samples': 5000, 'unlisted': 0, 'chromosomes': [
        BLOCK22 + ['--ids', 'dot'],
        ['--chrom', 'X', '--transcripts', '20', '--ids', 'rs', '--males', 'mirrored'],
    ]},
    'n10000': {'samples': 10000, 'unlisted': 0, 'chromosomes': [BLOCK22 + ['--ids', 'ambiguous']]},
    'whole22': {'samples': 2500, 'unlisted': 0, 'chromosomes': [
        ['--chrom', '21', '--transcripts', '20', '--ids', 'rs'],
        ['--chrom', '22', '--transcripts', '0', '--ids', 'dot'],
    ]},
    'dupseq': {'fasta': ['--records', '400000']},
    # large tier: 26,000 samples (32-row chunks); 100,000 samples on 3 transcripts per chromosome (8-row chunks);
    # TTN, the longest transcript (363 exons, 109 kb)
    'n26000': {'samples': 26000, 'unlisted': 0, 'chromosomes': [BLOCK22 + ['--ids', 'dot']]},
    # (three short transcripts with both strands per chromosome: legacy prohap.py holds a per-transcript table of
    # rows x 100,000 samples in memory, about 6 GB for the 240 rows of the largest here)
    'n100000': {'samples': 100000, 'unlisted': 0, 'chromosomes': [
        ['--chrom', '21', '--first', 'ENST00000290310', '--transcripts', '3', '--ids', 'ambiguous'],
        ['--chrom', '22', '--first', 'ENST00000290765', '--transcripts', '3', '--ids', 'dot'],
    ]},
    'ttn': {'samples': 2500, 'unlisted': 0, 'chromosomes': [
        ['--chrom', '2', '--first', 'ENST00000589042', '--transcripts', '1', '--ids', 'dot'],
    ]},
}


def default_refs():
    root = os.environ.get('PROHAP_REFS', '')
    found = sorted(d for d in (os.listdir(root) if root and os.path.isdir(root) else []) if d.startswith('ensembl.v'))
    return os.path.join(root, found[0]) if len(found) == 1 else None


def sha256(path):
    with open(path, 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def make(name, refs):
    spec = DATASETS[name]
    with open(os.path.join(refs, 'release.txt')) as fh:
        release = fh.read().strip()
    tool = SIMULATE_FASTA if 'fasta' in spec else SIMULATE
    stamp = {'dataset': spec, 'simulator_sha256': sha256(tool), 'ensembl_release': release}
    path = os.path.join(CACHE, name)
    try:
        with open(os.path.join(path, 'stamp.json')) as fh:
            if json.load(fh) == stamp:
                return path
    except (OSError, ValueError):
        pass
    print('Simulating fixture data %s into %s' % (name, path), file=sys.stderr, flush=True)
    partial = path + '.partial'
    shutil.rmtree(partial, ignore_errors=True)
    os.makedirs(partial)
    if 'fasta' in spec:
        commands = [[sys.executable, SIMULATE_FASTA, '--refs', refs] + spec['fasta']
                    + ['--out', os.path.join(partial, 'input.fa.gz')]]
    else:
        commands = [[sys.executable, SIMULATE, '--refs', refs, '--samples', str(spec['samples']), '--unlisted',
                     str(spec['unlisted']), '--seed', '0'] + args + ['--out', partial] for args in spec['chromosomes']]
    for command in commands:
        subprocess.run(command, check=True)
    with open(os.path.join(partial, 'stamp.json'), 'w') as fh:
        json.dump(stamp, fh, indent=1)
        fh.write('\n')
    shutil.rmtree(path, ignore_errors=True)
    os.rename(partial, path)
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('dataset', nargs='?', choices=sorted(DATASETS))
    parser.add_argument('--list', action='store_true', help='list the datasets')
    parser.add_argument('--refs', default=default_refs(),
                        help='references directory ensembl.v<release> (default: the one in $PROHAP_REFS)')
    args = parser.parse_args(argv)
    if args.list:
        print('\n'.join(sorted(DATASETS)))
        return
    if not args.dataset:
        parser.error('name a dataset (--list shows them)')
    if not args.refs:
        parser.error('--refs is required ($PROHAP_REFS does not hold exactly one ensembl.v<release>)')
    print(make(args.dataset, os.path.abspath(os.path.expanduser(args.refs))))


if __name__ == '__main__':
    main()
