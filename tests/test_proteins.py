"""src/modules/proteins.py (R2) against the reference: translate, reverse_complement and get_affected_codons against
Biopython and legacy coordinates_toolbox, and process_transcript against legacy process_haplotypes on the
reference's own haplotype table for crafted transcripts and variants, with every combination of settings.

    pytest tests/test_proteins.py
"""
import io
import random
import warnings

import pytest

pd = pytest.importorskip('pandas')

from Bio import BiopythonWarning  # noqa: E402
from Bio.Seq import Seq  # noqa: E402

from fast_cases import make_protein_case  # noqa: E402
from modules import proteins  # noqa: E402
from modules.legacy_import import legacy  # noqa: E402

warnings.simplefilter('ignore', BiopythonWarning)


class SerialPool:
    def __init__(self, processes=None):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def map(self, func, items):
        return [func(item) for item in items]


def test_translate_matches_biopython():
    rng = random.Random(1)
    for _ in range(3000):
        n = rng.choice([0, 1, 2, 3, 5, 30, 89, 90, 91, 92, 93, 400])
        alphabet = rng.choice(['ACGT', 'ACGT', 'ACGTN', 'acgt', 'ACGTRY'])
        seq = ''.join(rng.choice(alphabet) for _ in range(n))
        assert proteins.translate(seq) == str(Seq(seq).transcribe().translate()), seq


def test_reverse_complement_matches_biopython():
    rng = random.Random(2)
    for _ in range(500):
        seq = ''.join(rng.choice(rng.choice(['ACGT', 'ACGTN'])) for _ in range(rng.randint(0, 20)))
        assert proteins.reverse_complement(seq) == str(Seq(seq).reverse_complement())


def test_get_affected_codons_matches_legacy():
    rng = random.Random(3)
    cdna = ''.join(rng.choice('ACGT') for _ in range(200))
    for _ in range(5000):
        args = (rng.randint(-6, 206), rng.randint(0, 12), rng.choice([-1, 0, 1, 2]), rng.randint(0, 20))
        expected = legacy.coordinates_toolbox.get_affected_codons(Seq(cdna), *args)
        assert proteins.get_affected_codons(cdna, *args) == expected, args


def test_add_population_freqs_matches_legacy():
    rng = random.Random(4)
    for _ in range(3000):
        strings = []
        for _ in range(rng.randint(1, 6)):
            codes = rng.sample(['GBR', 'FIN', 'YRI', 'CHB', 'PEL', 'X1'], rng.randint(1, 5))
            strings.append(';'.join(c + ':{:.5f}'.format(rng.choice([rng.random(), rng.random() / 1000, 0.5, 0.0]))
                                    for c in codes))
        expected = strings[0]
        for s in strings[1:]:
            expected = legacy.process_haplotypes.add_population_freqs(expected, s)
        assert proteins.add_population_freqs(strings) == expected, strings


class TableRow:
    """A row of the reference's haplotype table, with the interface of haplotypes.Row."""
    def __init__(self, row):
        self.changes = row['Changes']
        self.afs = row['AlleleFrequencies']
        self.ids = row['VCF_IDs']
        self.count = row['Count']
        self.row = row

    def frequency(self, cohort):
        return self.row['Frequency']

    def frequency_strings(self, cohort):
        return [self.row['Frequency_population'], self.row['Frequency_superpopulation']]

    def sample_string(self, cohort):
        return self.row['Samples']


def reference_table(seed, tmp_path, monkeypatch):
    features, db, cdnas, samples, text = make_protein_case(seed)
    tables = tmp_path / 'tables'
    tables.mkdir(exist_ok=True)
    colnames = legacy.vcf_reader.parse_vcf(features, io.StringIO(text), db, 0, str(tables))
    sample_ids = [s for s in colnames if s.startswith('S')]
    metadata = pd.DataFrame({'Sample name': sample_ids, 'Sex': 'female',
                             'Population code': ['P%d' % (i % 3) for i in range(len(sample_ids))],
                             'Superpopulation code': ['Q%d' % (i % 2) for i in range(len(sample_ids))]})
    monkeypatch.setattr(legacy.get_haplotypes, 'Pool', SerialPool)
    haplo_df = legacy.get_haplotypes.get_gene_haplotypes(features, sample_ids, str(tables),
                                                         str(tmp_path / 'log'), 1, False, 1, 2, metadata)
    return features, db, cdnas, haplo_df


SETTINGS = [  # force_rf, min_freq, min_count, ignore_UTR, skip_start_lost, output cDNAs
    (1, -1, 0, 1, 1, True),
    (1, -1, 0, 0, 0, True),
    (1, 0.05, 0, 1, 0, False),
    (0, -1, 3, 0, 1, True),
    (0, -1, 0, 1, 0, False),
]


@pytest.mark.parametrize('seed', [1, 2, 3, 4])
@pytest.mark.parametrize('settings', SETTINGS, ids=['tre', 'utr-startlost', 'freq', 'count-norf', 'norf'])
def test_process_transcript_matches_process_haplotypes(seed, settings, tmp_path, monkeypatch):
    force_rf, min_freq, min_count, ignore_utr, skip_start_lost, output_cdnas = settings
    features, db, cdnas, haplo_df = reference_table(seed, tmp_path, monkeypatch)
    monkeypatch.setattr(legacy.process_haplotypes, 'Pool', SerialPool)
    expected, expected_proteins, expected_cdnas = legacy.process_haplotypes.process_haplotypes(
        features, haplo_df, cdnas, db, '22', 'haplo_chr22', force_rf, 1, min_freq, min_count, ignore_utr,
        skip_start_lost, output_cdnas)

    options = proteins.Settings('22', force_rf, min_freq, min_count, ignore_utr, skip_start_lost, output_cdnas)
    rows, protein_seqs, cdna_seqs = [], [], []
    name = lambda i: 'haplo_chr22_' + hex(i)[2:]       # noqa: E731
    for f in features:
        exons = db.children(f, 'exon', 'start')
        starts = db.children(f, 'start_codon', 'start')
        stops = db.children(f, 'stop_codon', 'start')
        biotype = f.biotype if f.biotype is not None else '-'
        t = proteins.Transcript(f.id, f.start, f.end, f.strand, exons, starts[0].start if starts else None,
                                stops[0].start if stops else None, biotype,
                                cdnas[f.id]['sequence'] if f.id in cdnas else None)
        group = haplo_df[haplo_df['TranscriptID'] == f.id]
        table, prots, cds = proteins.process_transcript(t, [(i, TableRow(r)) for i, r in group.iterrows()], options,
                                                        None)
        for fields in table:
            fields[3] = name(fields[3])
        rows.extend(table)
        protein_seqs.extend(prots)
        cdna_seqs.extend(cds)

    got = pd.DataFrame(rows, columns=legacy.process_haplotypes.result_columns)
    assert got.to_csv(sep='\t', index=False) == expected.to_csv(sep='\t', index=False)
    assert [([name(h) for h in s['haplotypes']], s['rfs'], s['start'], s['sequence']) for s in protein_seqs] == \
        [(s['haplotypes'], s['rfs'], s['start'], str(s['sequence'])) for s in expected_proteins]
    assert [([name(h) for h in s['haplotypes']], s['start'], s['sequence']) for s in cdna_seqs] == \
        [(s['haplotypes'], s['start'], str(s['sequence'])) for s in expected_cdnas]
    assert len(expected) > (0 if min_freq > 0 or min_count > 0 else 5)
