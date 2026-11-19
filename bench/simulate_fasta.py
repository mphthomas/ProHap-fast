#!/usr/bin/env python3
"""Simulate the protein FASTA that merge_duplicate_seq.py reads in the workflow, at any size.

    python3 bench/simulate_fasta.py --refs ~/prohap-refs/ensembl.v113 --records 400000 --out input.fa.gz

In the workflow, merge_duplicate_seq.py reads the cleaned reference proteome, the tagged contaminants and the cleaned
haplotype FASTA, concatenated; at cohort scale the haplotype records dominate. This writes that file (gzip level 1):

- the reference: the MANE Select proteins of the references directory, with the header remove_stop_codons.py gives a
  protein without stop codons (`generic_ensref|ENST...|position_within_protein:0 pep ...`); proteins with a `*` are
  left out;
- the contaminants: crap_tagged.fa.gz as it is;
- haplotype records up to --records in all, in prohap.py's format after remove_stop_codons.py
  (`generic_enshap|enshap_<chr>_<hex>|position_within_protein:0 matching_proteins:haplo_chr<chr>_<hex>;...
  start:0 reading_frame:0;...`): 85% a reference protein with 1-3 substitutions, 10% the sequence of an earlier
  haplotype record, 3% the reference sequence itself, and 2% short 3' UTR fragments (`..._3UTR_1`, 6-20 residues from
  a small pool), so merge_duplicate_seq.py has duplicates to merge.

The same arguments write the same bytes.
"""
import argparse
import gzip
import io
import os
import sys

import numpy as np

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


def default_refs():
    root = os.environ.get('PROHAP_REFS', '')
    found = sorted(d for d in (os.listdir(root) if root and os.path.isdir(root) else []) if d.startswith('ensembl.v'))
    return os.path.join(root, found[0]) if len(found) == 1 else None


def read_fasta(fh):
    header, parts = None, []
    for line in fh:
        if line.startswith('>'):
            if header is not None:
                yield header, ''.join(parts)
            header, parts = line[1:].rstrip('\n'), []
        else:
            parts.append(line.strip())
    if header is not None:
        yield header, ''.join(parts)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--refs', default=default_refs(),
                        help='references directory ensembl.v<release> (default: the one in $PROHAP_REFS)')
    parser.add_argument('--records', type=int, default=400000, help='records in all (default 400,000)')
    parser.add_argument('--seed', type=int, default=0, help='random seed (default 0)')
    parser.add_argument('--out', required=True, help='output FASTA (.gz)')
    args = parser.parse_args(argv)
    if not args.refs:
        parser.error('--refs is required ($PROHAP_REFS does not hold exactly one ensembl.v<release>)')
    refs = os.path.abspath(os.path.expanduser(args.refs))
    with open(os.path.join(refs, 'release.txt')) as fh:
        release = fh.read().strip()
    rng = np.random.default_rng(np.random.SeedSequence([args.seed, 3]))

    with open(os.path.join(refs, 'transcripts_reference_%s.csv' % release)) as fh:
        mane = {line.split(',')[1].strip() for line in list(fh)[1:]}
    reference = []          # (header, sequence, chromosome)
    with open(os.path.join(refs, 'ensembl_reference_proteinDB_%s_tagged.fa' % release)) as fh:
        for header, seq in read_fasta(fh):
            tag, accession, description = header.split('|', 2)
            if accession in mane and '*' not in seq:
                chrom = description.split('chromosome:', 1)[1].split(':')[1] if 'chromosome:' in description else '1'
                reference.append(('%s|%s|position_within_protein:0 %s' % (tag, accession, description), seq, chrom))
    with gzip.open(os.path.join(refs, 'crap_tagged.fa.gz'), 'rt') as fh:
        contaminants = list(read_fasta(fh))
    n_haplo = args.records - len(reference) - len(contaminants)
    if n_haplo < 0:
        parser.error('--records must be at least %d (reference and contaminants)'
                     % (len(reference) + len(contaminants)))

    counters = {}

    def next_id(chrom):
        counters[chrom] = counters.get(chrom, -1) + 1
        return '%x' % counters[chrom]

    fragments = [''.join(rng.choice(list(AMINO_ACIDS), size=int(rng.integers(6, 21)))) for _ in range(200)]
    kinds = rng.random(n_haplo)
    picks = rng.integers(0, len(reference), n_haplo)
    written = []            # sequences of the haplotype records so far
    with open(args.out, 'wb') as raw, gzip.GzipFile(filename='', mode='wb', compresslevel=1, fileobj=raw,
                                                    mtime=0) as gz, io.TextIOWrapper(gz, newline='\n') as out:
        for header, seq, _ in reference:
            out.write('>%s\n%s\n' % (header, seq))
        for header, seq in contaminants:
            out.write('>%s\n%s\n' % (header, seq))
        for i in range(n_haplo):
            _, ref_seq, chrom = reference[picks[i]]
            accession = 'enshap_%s_%s' % (chrom, next_id(chrom))
            position = 0
            if kinds[i] < 0.85:
                seq = list(ref_seq)
                for _ in range(int(rng.integers(1, 4))):
                    k = int(rng.integers(0, len(seq)))
                    seq[k] = str(rng.choice([a for a in AMINO_ACIDS if a != seq[k]]))
                seq = ''.join(seq)
            elif kinds[i] < 0.95 and written:
                seq = written[int(rng.integers(0, len(written)))]
            elif kinds[i] < 0.98:
                seq = ref_seq
            else:
                seq = fragments[int(rng.integers(0, len(fragments)))]
                accession += '_3UTR_1'
                position = int(rng.integers(100, 1000))
            written.append(seq)
            haplotypes = ['haplo_chr%s_%s' % (chrom, '%x' % int(rng.integers(0, 1 << 20)))
                          for _ in range(int(rng.integers(1, 4)))]
            out.write('>generic_enshap|%s|position_within_protein:%d matching_proteins:%s start:0 '
                      'reading_frame:%s\n%s\n'
                      % (accession, position, ';'.join(haplotypes), ';'.join('0' * len(haplotypes)), seq))
    print('[simulate_fasta] %d records: %d reference, %d contaminants, %d haplotype -> %s'
          % (args.records, len(reference), len(contaminants), n_haplo, args.out), file=sys.stderr)


if __name__ == '__main__':
    main()
