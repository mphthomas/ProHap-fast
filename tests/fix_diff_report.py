#!/usr/bin/env python3
"""What each fix of prohap.py -fix_bugs 1 changes in the output, fix by fix.

    python3 tests/fix_diff_report.py [--work DIR] [--keep] [--chromosomes 21,22] [DATASET...]   # default 1kg

For each chromosome of the dataset (tests/fix_oracle.py has the datasets) it runs src/prohap.py with no fix, then
with the fixes added one at a time in the order of modules/fixes.ALL (PROHAP_FIXES; with none it writes the
reference's output, test_fix_bugs.py checks that), every haplotype kept, and compares each run with the one before.
Each difference is so attributed to the one fix that was added. Per fix and in total:

    removed copies       haplotype copies listed in prohap_chr*.log: no longer removed / newly removed
    rows                 haplotype rows, identified by transcript, DNA_changes and VCF_IDs: gained / lost /
                         same row with other samples
    FASTA records        distinct protein sequences per transcript: gained / lost
    copies               haplotype copies by outcome before -> after: removed, ref (no row), a protein
"""
import argparse
import collections
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))
import fix_oracle  # noqa: E402
from modules import fixes  # noqa: E402


def read_run(out):
    """(removed copies, rows, FASTA records, outcome per copy) of one prohap.py run."""
    out = Path(out)
    sequences = {}
    lines = (out / 'haplo.fa').read_text().split('\n')
    for header, seq in zip(lines[0::2], lines[1::2]):
        for h in header.split('matching_proteins:')[1].split()[0].split(';'):
            sequences[h] = seq
    rows = {}
    outcome = {}
    records = set()
    with open(out / 'haplo.tsv') as fh:
        columns = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            row = dict(zip(columns, line.rstrip('\n').split('\t')))
            key = (row['TranscriptID'], row['DNA_changes'], row['VCF_IDs'])
            rows[key] = row['samples']
            seq = sequences.get(row['HaplotypeID'])
            if seq is not None:
                records.add((row['TranscriptID'], seq))
            for copy in row['samples'].split(';'):
                outcome[(row['TranscriptID'], copy)] = seq
    removed = set()
    with open(out / 'prohap.log') as fh:
        fh.readline()
        for line in fh:
            transcript, _, copies = line.rstrip('\n').split('\t')
            for copy in (copies.split(';') if copies else []):
                removed.add((transcript, copy))
                outcome[(transcript, copy)] = 'removed'
    return removed, rows, records, outcome


def kind(x):
    return 'ref' if x is None else ('removed' if x == 'removed' else 'protein')


def diff(before, after):
    removed_a, rows_a, records_a, outcome_a = before
    removed_b, rows_b, records_b, outcome_b = after
    counts = collections.Counter()
    counts['removed copies: no longer removed'] = len(removed_a - removed_b)
    counts['removed copies: newly removed'] = len(removed_b - removed_a)
    counts['rows: gained'] = len(set(rows_b) - set(rows_a))
    counts['rows: lost'] = len(set(rows_a) - set(rows_b))
    counts['rows: same row, other samples'] = sum(1 for k in set(rows_a) & set(rows_b) if rows_a[k] != rows_b[k])
    counts['FASTA records: gained'] = len(records_b - records_a)
    counts['FASTA records: lost'] = len(records_a - records_b)
    for copy in set(outcome_a) | set(outcome_b):
        a, b = outcome_a.get(copy), outcome_b.get(copy)
        if a != b:
            ka, kb = kind(a), kind(b)
            counts['copies: %s -> %s' % (ka, 'other protein' if ka == kb == 'protein' else kb)] += 1
    return counts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('datasets', nargs='*', default=['1kg'])
    parser.add_argument('--work')
    parser.add_argument('--keep', action='store_true')
    parser.add_argument('--chromosomes', help='only these chromosomes, comma-separated')
    args = parser.parse_args(argv)
    work = Path(args.work or tempfile.mkdtemp(prefix='fix_diff_report.')).resolve()
    work.mkdir(parents=True, exist_ok=True)
    known = fix_oracle.datasets(work)
    levels = [tuple(fixes.ALL[:k]) for k in range(len(fixes.ALL) + 1)]
    names = ['+' + level[-1] if level else 'none' for level in levels]
    totals = {name: collections.Counter() for name in names[1:]}
    baseline = collections.Counter()
    for dataset in args.datasets:
        for chrom, vcf_glob, af, samples, transcripts in known[dataset]():
            if args.chromosomes and chrom not in args.chromosomes.split(','):
                continue
            runs = []
            for name, level in zip(names, levels):
                out = work / dataset / ('chr' + chrom) / name
                fix_oracle.run_prohap(out, chrom, vcf_glob, af, samples, transcripts, fixes=level)
                runs.append(read_run(out))
            removed, rows, records, _ = runs[0]
            baseline.update({'removed copies': len(removed), 'rows': len(rows), 'FASTA records': len(records)})
            for k in range(1, len(runs)):
                counts = diff(runs[k - 1], runs[k])
                totals[names[k]].update(counts)
                print('%s chr%s, %s: %s' % (dataset, chrom, names[k],
                                            ', '.join('%s %d' % kv for kv in sorted(counts.items()) if kv[1])),
                      flush=True)
    print()
    print('Without fixes: %s' % ', '.join('%s %d' % kv for kv in sorted(baseline.items())))
    keys = sorted(set().union(*totals.values()))
    width = max(len(k) for k in keys) if keys else 10
    print('%-*s  %s' % (width, '', '  '.join('%12s' % n for n in names[1:])))
    for key in keys:
        print('%-*s  %s' % (width, key, '  '.join('%12d' % totals[n][key] for n in names[1:])))
    if not args.keep and not args.work:
        shutil.rmtree(work)


if __name__ == '__main__':
    main()
