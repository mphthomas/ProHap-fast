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
