"""Write transcripts.csv, the transcript list of the 1kG fixture: a small part of chromosomes 1, 21, 22 and X that
still exercises everything, so a legacy run takes seconds instead of minutes.

The legacy code spends about 0.1 s per transcript whether or not it has variants (it writes and re-reads a table with
one column per sample for every transcript), and only about one MANE transcript in five overlaps a variant of the 1kG
sample data. So the list keeps, on each VCF chunk (chromosome 1 comes in two), 20 transcripts whose exons hold at
least one record and 5 whose exons hold none, evenly spread along the chromosome. Chromosome X has only 7 records,
all in PAR1, and 4 transcripts that overlap them. Rows keep the order of the reference list.

    PYTHONHASHSEED=0 python3 tests/fixtures/1kg/make_transcripts.py ~/prohap-refs/ensembl.v113 \\
        > tests/fixtures/1kg/transcripts.csv

Rerun it for a new Ensembl release; the outputs of the fixture change with it.
"""
import bisect
import os
import sys
import tarfile

import gffutils

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
TARBALL = os.path.join(ROOT, 'sample_data', 'sample_1kGP_common_global.tar.gz')
CHUNKS = {'1': ['sample_chr1_1kGP_GRCh38_phased.vcf', 'sample_chr1_1kGP_GRCh38_phased_2.vcf'],
          '21': ['sample_chr21_1kGP_GRCh38_phased.vcf'],
          '22': ['sample_chr22_1kGP_GRCh38_phased.vcf'],
          'X': ['sample_chrX_1kGP_GRCh38_phased.vcf']}
WITH, WITHOUT = 20, 5


def positions(tar, name):
    """POS of every record of one VCF chunk, sorted."""
    with tar.extractfile(name) as f:
        return sorted(int(line.split(b'\t', 2)[1]) for line in f if not line.startswith(b'#'))


def spread(items, n):
    """n items evenly spaced through items (all of them if there are no more than n)."""
    if len(items) <= n:
        return items
    return [items[(2 * i + 1) * len(items) // (2 * n)] for i in range(n)]


def main(refs):
    release = open(os.path.join(refs, 'release.txt')).read().strip()
    with open(os.path.join(refs, 'transcripts_reference_%s.csv' % release)) as f:
        reference = [line.rstrip('\n').split(',') for line in f][1:]
    keep = set()
    with tarfile.open(TARBALL) as tar:
        for chrom, chunks in CHUNKS.items():
            db = gffutils.FeatureDB(os.path.join(refs, 'Homo_sapiens.GRCh38.%s.chr_patch_hapl_scaff_chr%s.db'
                                                 % (release, chrom)))
            transcripts = sorted((db[t].start, t) for c, t in reference if c == chrom)
            for chunk in chunks:
                pos = positions(tar, chunk)
                # with two chunks, each one takes the transcripts that start within its records
                inside = [t for start, t in transcripts if len(chunks) == 1 or pos[0] <= start <= pos[-1]]

                def overlaps(t):
                    for exon in db.children(t, featuretype='exon'):
                        i = bisect.bisect_left(pos, exon.start)
                        if i < len(pos) and pos[i] <= exon.end:
                            return True
                    return False

                hit = [t for t in inside if overlaps(t)]
                miss = [t for t in inside if t not in hit]
                keep.update(spread(hit, WITH), spread(miss, WITHOUT))
    print('chromosome,transcriptID')
    for chrom, t in reference:
        if t in keep:
            print(chrom + ',' + t)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit('usage: make_transcripts.py REFS_DIR (e.g. ~/prohap-refs/ensembl.v113)')
    main(os.path.expanduser(sys.argv[1]))
