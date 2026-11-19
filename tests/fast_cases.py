"""Crafted inputs shared by the tests of the fast prohap.py modules (test_sweep, test_haplotypes, test_proteins):
transcripts with exons, a VCF, and the reference's parse_vcf tables for them."""
import gzip
import io
import random
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT / 'src') not in sys.path:
    sys.path.insert(0, str(ROOT / 'src'))
from modules import sweep  # noqa: E402
from modules.legacy_import import legacy  # noqa: E402

COLUMNS = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']


class FakeDB:
    def __init__(self, exons):
        self.exons = exons

    def children(self, transcript, featuretype, order_by):
        assert (featuretype, order_by) == ('exon', 'start')
        return [SimpleNamespace(start=s, end=e) for s, e in self.exons[transcript.id]]


def make_case(seed, n_samples=12, n_records=300, unsorted=False):
    rng = random.Random(seed)
    exons = {}
    for i in range(20):
        start = rng.randint(1, 40) * 50
        if i % 5 == 0 and exons:                     # same start as an earlier transcript
            start = exons[rng.choice(sorted(exons))][0][0]
        blocks = []
        pos = start
        for _ in range(rng.randint(1, 4)):
            length = rng.randint(5, 150)
            blocks.append((pos, pos + length))
            pos += length + rng.randint(1, 40)
        exons['T%02d' % i] = blocks
    exons['TNOEX'] = []
    transcripts = [SimpleNamespace(id=t, start=e[0][0] if e else 900, end=e[-1][1] if e else 950)
                   for t, e in exons.items()]
    transcripts.sort(key=lambda t: t.start)
    samples = ['S%02d' % i for i in range(n_samples)]
    lines = ['##fileformat=VCFv4.2', '##contig=<ID=22>', '\t'.join(COLUMNS + samples)]
    positions = sorted(rng.randint(1, 2600) for _ in range(n_records))
    if unsorted:                                    # a few records out of order
        for _ in range(8):
            i = rng.randrange(n_records - 3)
            j = i + rng.randint(1, 3)
            positions[i], positions[j] = positions[j], positions[i]
    for n, pos in enumerate(positions):
        ref, alt = rng.choice([('A', 'G'), ('C', 'T'), ('ACGTACGTAC', 'A'), ('G', 'GTT'), ('N', 'A'), ('A', '*'),
                               ('T', '-'), ('-', 'A')])
        vcf_id = rng.choice(['.', '.', '.', 'rs%d' % n, '2e5', 'NA', 'True', '0012', '7'])
        info = rng.choice(['AC=1;AN=24;AF=0.1', 'MAF=0.3', '.', 'DP=4'])
        tokens = []
        for _ in samples:
            a, b = rng.choice('0000001.'), rng.choice('0000001.')
            tokens.append(a + '|' + b)
        lines.append('\t'.join(['22', str(pos), vcf_id, ref, alt, '.', 'PASS', info, 'GT'] + tokens))
    return transcripts, exons, samples, '\n'.join(lines) + '\n'


def reference(transcripts, exons, text, tmp_path, monkeypatch):
    """parse_vcf's tables, read back as get_haplotypes reads them, in the order parse_vcf wrote them."""
    order = []
    original = pd.DataFrame.to_csv

    def to_csv(self, path, *args, **kwargs):
        order.append(Path(path).stem)
        return original(self, path, *args, **kwargs)
    monkeypatch.setattr(pd.DataFrame, 'to_csv', to_csv)
    out = tmp_path / 'legacy_tables'
    out.mkdir()
    colnames = legacy.vcf_reader.parse_vcf(transcripts, io.StringIO(text), FakeDB(exons), 0, str(out))
    monkeypatch.setattr(pd.DataFrame, 'to_csv', original)
    tables = {}
    for t in transcripts:
        df = pd.read_csv(out / (t.id + '.tsv'), sep='\t')
        tables[t.id] = df
    return colnames, order, tables


def fast(transcripts, exons, samples, included, text, tmp_path):
    path = tmp_path / 'in.vcf.gz'
    with gzip.open(path, 'wt') as f:
        f.write(text)
    reader = sweep.VCFReader(str(path))
    first_number, columns, first = sweep.read_header(reader)
    rank = [included.index(s) if s in included else -1 for s in columns[9:]]
    parser = sweep.RecordParser(len(columns) - 9, rank)
    infos = []
    for t in transcripts:
        starts, max_ends = sweep.exon_index(exons[t.id])
        infos.append(SimpleNamespace(start=t.start, end=t.end, exon_starts=starts, exon_max_ends=max_ends))
    released = list(sweep.assign_records(reader, first, first_number, infos, parser))
    reader.finish()
    return columns, released


# ---------------------------------------------------------------------------------------------------------------
# a genome, transcripts with their cDNA, start and stop codons, and variants whose REF alleles match the genome

_COMPLEMENT = str.maketrans('ACGT', 'TGCA')


class Feature(SimpleNamespace):
    """A gffutils transcript feature as process_haplotypes uses it: attributes and ['transcript_biotype']."""
    def __getitem__(self, key):
        if key == 'transcript_biotype' and self.biotype is not None:
            return [self.biotype]
        raise KeyError(key)


class FullDB:
    def __init__(self, children):
        self.data = children          # {transcript id: {featuretype: [(start, end)]}}

    def children(self, transcript, featuretype, order_by):
        assert order_by == 'start'
        return [SimpleNamespace(start=s, end=e) for s, e in sorted(self.data[transcript.id].get(featuretype, []))]


def make_protein_case(seed, n_samples=40, n_transcripts=8):
    """(features, db, cdnas, samples, vcf text): transcripts on both strands with and without start and stop codons,
    SNVs, deletions (some across exon ends, some with an SNV inside on the same copies), insertions, multi-allelic
    sites written as split records, a few REF alleles that do not match the genome and '-' alleles, and genotypes
    drawn from a few founder haplotypes so that copies share haplotypes."""
    rng = random.Random(seed)
    genome = ''.join(rng.choice('ACGT') for _ in range(6000))
    base = lambda pos: genome[pos - 1]                      # noqa: E731  (1-based)
    features = []
    children = {}
    cdnas = {}
    position = 100
    for i in range(n_transcripts):
        exons = []
        start = position
        for _ in range(rng.randint(1, 5)):
            length = rng.randint(20, 160)
            exons.append((position, position + length - 1))
            position += length + rng.randint(5, 60)
        if rng.random() < 0.3:
            position = start + rng.randint(10, 200)        # overlapping transcripts
        tid = 'T%02d' % i
        strand = rng.choice('+-')
        spliced = ''.join(genome[s - 1:e] for s, e in exons)
        cdna = spliced if strand == '+' else spliced.translate(_COMPLEMENT)[::-1]
        if rng.random() < 0.15:
            cdna = None                                     # not in the cDNA FASTA
        codons = {'exon': exons}
        coding = [p for s, e in exons for p in range(s, e + 1)]
        mode = rng.choice(['start', 'start', 'start', 'stop', 'none'])
        if mode == 'start':
            codons['start_codon'] = [(c, c + 2) for c in [rng.choice(coding[:len(coding) // 2])]]
            codons['stop_codon'] = [(c, c + 2) for c in [rng.choice(coding[len(coding) // 2:])]]
        elif mode == 'stop':
            codons['stop_codon'] = [(c, c + 2) for c in [rng.choice(coding)]]
        children[tid] = codons
        if cdna is not None:
            cdnas[tid] = {'sequence': cdna}
        features.append(Feature(id=tid, start=exons[0][0], end=exons[-1][1], strand=strand,
                                biotype=rng.choice(['protein_coding', None])))
    features.sort(key=lambda f: f.start)

    samples = ['S%03d' % i for i in range(n_samples)]
    founders = [[rng.random() < 0.5 for _ in range(400)] for _ in range(6)]
    copy_founder = [rng.randrange(6) for _ in range(2 * n_samples)]
    lines = ['##fileformat=VCFv4.2', '\t'.join(COLUMNS + samples)]
    covered = sorted({p for t in children.values() for s, e in t['exon'] for p in range(s - 3, e + 4)})
    sites = sorted(rng.sample(covered, min(len(covered), 160)))
    for n, pos in enumerate(sites):
        kind = rng.random()
        ref = base(pos)
        if kind < 0.04:                                 # a REF that does not match the genome (nor the cDNA)
            ref = rng.choice([b for b in 'ACGT' if b != ref])
            alts = [rng.choice([b for b in 'ACGT' if b != ref])]
        elif kind < 0.07:                               # '-' alleles
            ref, alts = rng.choice([(ref, ['-']), ('-', [rng.choice('ACGT')])])
        elif kind < 0.6:
            alts = [rng.choice([b for b in 'ACGT' if b != ref])]
        elif kind < 0.75:
            length = rng.randint(1, 8)
            ref = genome[pos - 1:pos + length]
            alts = [ref[0]]
        elif kind < 0.9:
            alts = [ref + ''.join(rng.choice('ACGT') for _ in range(rng.randint(1, 6)))]
        else:
            alts = rng.sample([b for b in 'ACGT' if b != ref], 2)
        common = rng.random() < 0.4
        for j, alt in enumerate(alts):
            tokens = []
            for s in range(n_samples):
                pair = []
                for h in range(2):
                    c = 2 * s + h
                    carries = founders[copy_founder[c]][n % 400] if common else rng.random() < 0.08
                    pair.append('1' if carries and (len(alts) == 1 or c % 2 == j) else '0')
                tokens.append('|'.join(pair))
            vcf_id = rng.choice(['.', '.', 'rs%d' % (n // 2), 'NA'])
            lines.append('\t'.join(['22', str(pos), vcf_id, ref, alt, '.', 'PASS',
                                    'AC=1;AN=%d;AF=%.3f' % (2 * n_samples, rng.random()), 'GT'] + tokens))
            if len(ref) > 1 and rng.random() < 0.3:          # an SNV inside the deletion on the same copies
                inner = pos + rng.randint(1, len(ref) - 1)
                lines.append('\t'.join(['22', str(inner), '.', base(inner), rng.choice(
                    [b for b in 'ACGT' if b != base(inner)]), '.', 'PASS', 'AF=0.1', 'GT'] + tokens))
    # keep the VCF sorted by position (the inner SNVs come after their deletion)
    body = sorted(lines[2:], key=lambda line: int(line.split('\t')[1]))
    text = '\n'.join(lines[:2] + body) + '\n'
    return features, FullDB(children), cdnas, samples, text
