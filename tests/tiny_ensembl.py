"""A tiny synthetic Ensembl release in the format of the four FTP downloads, for tests that need references.

Chromosomes 21 and 22 plus one scaffold. Chromosome 21 has a gene with a MANE Select transcript and a second
transcript; chromosome 22 and the scaffold each have a gene whose only tagged transcript is Ensembl canonical. So the
MANE list built from it is ENST00000000211 (21), ENST00000000221 (22) and ENST00000000991 (scaffold).
"""
import gzip
import os

RELEASE = '113'
ANNOTATION = 'Homo_sapiens.GRCh38.%s.chr_patch_hapl_scaff' % RELEASE
CHROMOSOMES = ['21', '22']
MANE_LIST = 'chromosome,transcriptID\n21,ENST00000000211\n22,ENST00000000221\nKI270728.1,ENST00000000991\n'

# (chromosome, gene, [(transcript, tags, [(exon start, exon end)])])
GENES = [
    ('21', 'ENSG00000000021', [('ENST00000000211', ['MANE_Select', 'Ensembl_canonical'], [(1000, 1200), (2000, 2300)]),
                               ('ENST00000000212', [], [(1000, 1200), (3000, 3100)])]),
    ('22', 'ENSG00000000022', [('ENST00000000221', ['Ensembl_canonical'], [(5000, 5300)])]),
    ('KI270728.1', 'ENSG00000000099', [('ENST00000000991', ['Ensembl_canonical'], [(100, 400)])]),
]


def _gtf_lines():
    yield '#!genome-build GRCh38.p14\n'
    yield '#!genome-version GRCh38\n'
    for chrom, gene, transcripts in GENES:
        start = min(s for _, _, exons in transcripts for s, _ in exons)
        end = max(e for _, _, exons in transcripts for _, e in exons)
        gattr = 'gene_id "%s"; gene_version "1"; gene_source "ensembl_havana"; gene_biotype "protein_coding";' % gene
        yield '\t'.join([chrom, 'ensembl_havana', 'gene', str(start), str(end), '.', '+', '.', gattr]) + '\n'
        for tr, tags, exons in transcripts:
            tattr = (gattr.rstrip(';') + '; transcript_id "%s"; transcript_version "1"; '
                     'transcript_biotype "protein_coding";' % tr)
            tattr += ''.join(' tag "%s";' % t for t in tags)
            yield '\t'.join([chrom, 'ensembl_havana', 'transcript', str(exons[0][0]), str(exons[-1][1]),
                             '.', '+', '.', tattr]) + '\n'
            for n, (s, e) in enumerate(exons, 1):
                eattr = tattr + ' exon_number "%d"; exon_id "%s";' % (n, tr.replace('ENST', 'ENSE') + str(n))
                yield '\t'.join([chrom, 'ensembl_havana', 'exon', str(s), str(e), '.', '+', '.', eattr]) + '\n'


def _fasta(kind):
    for chrom, gene, transcripts in GENES:
        for tr, _, exons in transcripts:
            loc = 'chromosome:GRCh38:%s:%d:%d:1' % (chrom, exons[0][0], exons[-1][1])
            if kind == 'pep':
                yield '>%s.1 pep %s gene:%s.1 transcript:%s.1 gene_biotype:protein_coding ' \
                      'transcript_biotype:protein_coding gene_symbol:G%s\n' % (tr.replace('ENST', 'ENSP'), loc, gene,
                                                                               tr, tr[-3:])
                yield 'MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\nKKKKKKKKKKKKKKKK*\n'
            elif kind == 'cdna':
                yield '>%s.1 cdna %s gene:%s.1 gene_biotype:protein_coding\n' % (tr, loc, gene)
                yield 'ATG' + 'GCA' * 30 + 'TAA\n'
    if kind == 'ncrna':
        yield '>ENST00000000500.1 ncrna chromosome:GRCh38:21:9000:9100:1 gene:ENSG00000000500.1\nACGUACGU\n'


def write_downloads(dest):
    """Write the four gzipped downloads into dest; return dest."""
    os.makedirs(dest, exist_ok=True)
    files = {
        ANNOTATION + '.gtf.gz': _gtf_lines(),
        'Homo_sapiens.GRCh38.cdna.all.fa.gz': _fasta('cdna'),
        'Homo_sapiens.GRCh38.ncrna.fa.gz': _fasta('ncrna'),
        'Homo_sapiens.GRCh38.pep.all.fa.gz': _fasta('pep'),
    }
    for name, lines in files.items():
        with gzip.open(os.path.join(dest, name), 'wt') as fh:
            fh.writelines(lines)
    return dest
