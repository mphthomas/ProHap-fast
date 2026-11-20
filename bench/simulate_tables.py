#!/usr/bin/env python3
"""Simulate the per-chromosome haplotype tables that merge_tables.py and haplo_extract_sample_names.py read, at any
size, without running the workflow.

    python3 bench/simulate_tables.py --chromosomes 1,2,3,4 --out DIR [--samples 500000] [--transcripts 60]
                                     [--rows 1950] [--listed 370000]

Writes DIR/haplo_chr<c>.tsv.gz (gzip level 1) in prohap.py's format: the 21 columns of its table, --transcripts
transcripts of --rows haplotypes each per chromosome. The defaults are the 500,000-sample benchmark block of prohap.py
(IMPLEMENTATION_NOTES section 5: 117k rows and 293 MB of text, 83% of it `samples`, per 1/315 of the exome), so one
chromosome here is one block. Each transcript lists about --listed haplotype copies (`S0123456:1`) over its rows, a
few common haplotypes carrying most of them and most rows carrying one to a few. The other columns hold what
prohap.py writes, with the value kinds that make pandas' typing matter: frequencies of up to 17 digits, booleans,
numbers, and IDs that are numbers in some rows (`0012`, `2e5`).

The same arguments write the same bytes.
"""
import argparse
import gzip
import os
import sys

import numpy as np

COLUMNS = ['TranscriptID', 'chromosome', 'transcript_biotype', 'HaplotypeID', 'VCF_IDs', 'DNA_changes',
           'allele_frequencies', 'cDNA_changes', 'all_protein_changes', 'variant_types', 'protein_changes',
           'reading_frame', 'protein_prefix_length', 'start_missing', 'start_lost', 'splice_sites_affected',
           'occurrence_count', 'frequency', 'frequency_population', 'frequency_superpopulation', 'samples']
TOKEN = 10              # 'S0123456:1'
POPULATIONS = ['GBR', 'YRI', 'CHB', 'PEL', 'GIH', 'FIN', 'ESN', 'JPT']


def copy_counts(rng, rows, listed):
    """Copies per haplotype: heavy-tailed, at least one each, about `listed` in all."""
    weights = rng.pareto(0.9, rows) + 1.0
    counts = np.maximum(1, np.floor(weights / weights.sum() * listed)).astype(np.int64)
    return counts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--chromosomes', default='1', help='comma-separated (default 1)')
    parser.add_argument('--samples', type=int, default=500000, help='cohort size (default 500,000)')
    parser.add_argument('--transcripts', type=int, default=60, help='transcripts per chromosome (default 60)')
    parser.add_argument('--rows', type=int, default=1950, help='haplotypes per transcript (default 1,950)')
    parser.add_argument('--listed', type=int, default=370000,
                        help='haplotype copies listed per transcript (default 370,000)')
    parser.add_argument('--seed', type=int, default=0, help='random seed (default 0)')
    parser.add_argument('--out', required=True, help='output directory')
    args = parser.parse_args(argv)
    if args.samples > 10 ** 7:
        parser.error('--samples must be below 10,000,000 (sample names have 7 digits)')
    os.makedirs(args.out, exist_ok=True)

    for k, chrom in enumerate(args.chromosomes.split(',')):
        rng = np.random.default_rng(np.random.SeedSequence([args.seed, 5, k]))
        # one long run of random tokens; a row's samples are a slice of it (fixed-width tokens)
        pool_tokens = min(4 * args.listed, 2 * args.samples)
        ids = rng.integers(0, args.samples, pool_tokens)
        copies = rng.integers(1, 3, pool_tokens)
        pool = ';'.join('S%07d:%d' % (i, c) for i, c in zip(ids.tolist(), copies.tolist()))
        path = os.path.join(args.out, 'haplo_chr%s.tsv.gz' % chrom)
        n_rows = 0
        with open(path, 'wb') as raw, gzip.GzipFile(filename='', mode='wb', compresslevel=1, fileobj=raw,
                                                    mtime=0) as out:
            out.write(('\t'.join(COLUMNS) + '\n').encode())
            hap = 0
            for t in range(args.transcripts):
                transcript = 'ENST%011d' % (10 ** 6 * (k + 1) + t)
                counts = copy_counts(rng, args.rows, args.listed)
                hap += 1        # the reference haplotype's row, dropped
                lines = []
                n_variants = rng.integers(1, 4, args.rows)
                for r in range(args.rows):
                    n = int(n_variants[r])
                    count = int(min(counts[r], pool_tokens))
                    start = int(rng.integers(0, pool_tokens - count + 1))
                    samples = pool[start * (TOKEN + 1):(start + count) * (TOKEN + 1) - 1]
                    positions = rng.integers(10 ** 6, 10 ** 8, n)
                    afs = rng.random(n) / 3
                    vcf_ids = ';'.join(rng.choice(['rs%d' % rng.integers(10 ** 8), '%x' % rng.integers(1, 1 << 20),
                                                   '0012', '2e5']) for _ in range(n))
                    change = '%d:V>%d:D' % (r + 1, r + 1)
                    pops = ';'.join('%s:%.5f' % (p, rng.random() / 5)
                                    for p in POPULATIONS[:int(rng.integers(1, len(POPULATIONS) + 1))])
                    lines.append('\t'.join([
                        transcript, chrom, 'protein_coding', 'haplo_chr%s_%x' % (chrom, hap), vcf_ids,
                        ';'.join('%d:A>G' % p for p in positions.tolist()),
                        ';'.join(repr(round(a, int(rng.integers(2, 6)))) for a in afs.tolist()),
                        ';'.join('%d:A>G' % (p % 5000) for p in positions.tolist()), change,
                        'SAV' if r % 5 else 'synonymous', change if r % 5 else '', '0', str(int(rng.integers(0, 40))),
                        'False', 'True' if r % 97 == 0 else 'False', '-' if r % 13 else '3', str(count),
                        repr(count / (2 * args.samples)), pops, 'EUR:%r' % (rng.random() / 3), samples]))
                    hap += 1
                out.write(('\n'.join(lines) + '\n').encode())
                n_rows += len(lines)
        print('[simulate_tables] chr%s: %d rows -> %s' % (chrom, n_rows, path), file=sys.stderr)


if __name__ == '__main__':
    main()
