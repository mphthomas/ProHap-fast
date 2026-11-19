"""src/modules/haplotypes.py (R1b) against the reference's get_gene_haplotypes, on the reference's own per-transcript
tables for crafted VCFs: every row (changes, allele frequencies, VCF IDs, count, samples, frequency, population and
superpopulation frequencies) and the removed-samples log.

    pytest tests/test_haplotypes.py
"""
import random

import pytest

pd = pytest.importorskip('pandas')

from fast_cases import fast, make_case, reference  # noqa: E402
from modules import haplotypes  # noqa: E402
from modules.fallback import Unsupported  # noqa: E402
from modules.legacy_import import legacy  # noqa: E402
from modules.pandas_emulation import pandas_chunk_rows  # noqa: E402
from modules.sweep import transcript_columns  # noqa: E402

PAR1_TO, PAR2_FROM = 1000, 2000          # the crafted transcripts lie between 50 and 2,600


class SerialPool:
    def __init__(self, processes=None):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def map(self, func, items):
        return [func(item) for item in items]


def metadata(samples, rng):
    pops = {'GBR': 'EUR', 'FIN': 'EUR', 'YRI': 'AFR', 'CHB': 'EAS', 'PEL': 'AMR'}
    names = samples[:-2] + ['M01', 'M02']                       # two VCF samples missing, two extra rows
    return pd.DataFrame({
        'Sample name': names,
        'Sex': [rng.choice(['male', 'female', 'male']) for _ in names],
        'Population code': [rng.choice(sorted(pops)) for _ in names],
        'Superpopulation code': '',
    }).assign(**{'Superpopulation code': lambda d: d['Population code'].map(pops)})


def reference_haplotypes(transcripts, tmp_path, monkeypatch, samples_df, colnames, is_x):
    log = tmp_path / 'legacy.log'
    sample_ids = [s for s in colnames if s in set(samples_df['Sample name'])]
    monkeypatch.setattr(legacy.get_haplotypes, 'Pool', SerialPool)
    df = legacy.get_haplotypes.get_gene_haplotypes(transcripts, sample_ids, str(tmp_path / 'legacy_tables'),
                                                   str(log), 1, is_x, PAR1_TO, PAR2_FROM, samples_df.copy())
    return df, log.read_text(), sample_ids


def fast_haplotypes(transcripts, released, cohort, n_cols):
    rows = []
    log = ['TranscriptID\t#Removed\tRemovedSamples\n']
    chunk_rows = pandas_chunk_rows(n_cols)
    for index in sorted(i for i, _ in released):           # the reference's table follows all_transcripts
        records = dict(released)[index]
        t = transcripts[index]
        autosomal = cohort.is_autosomal(t.start, t.end)
        if records:
            ids = transcript_columns(records, chunk_rows)['ID']
            result, removed = haplotypes.transcript_rows(records, ids, autosomal, cohort)
        else:
            result, removed = [haplotypes.reference_only_row(cohort)], []
        for row in result:
            if row.members is None:
                # the reference haplotype: never written, so its samples are not computed
                rows.append([t.id, 'REF', '', '', row.count, None, row.frequency(cohort), None, None])
            else:
                rows.append([t.id, row.changes, row.afs, row.ids, row.count, row.sample_string(cohort),
                             row.frequency(cohort)] + row.frequency_strings(cohort))
        if result:
            labels = [cohort.labels[c] for c in list(removed)]
            log.append('%s\t%d\t%s\n' % (t.id, len(labels), ';'.join(labels)))
    return rows, ''.join(log)


@pytest.mark.parametrize('seed, is_x, n_samples', [(1, False, 30), (2, True, 30), (3, False, 30), (4, True, 30),
                                                   (5, False, 30), (6, True, 300), (7, False, 450)])
def test_rows_match_get_gene_haplotypes(seed, is_x, n_samples, tmp_path, monkeypatch):
    rng = random.Random(seed)
    transcripts, exons, samples, text = make_case(seed, n_samples=n_samples)
    samples_df = metadata(samples, rng)
    colnames, _, _ = reference(transcripts, exons, text, tmp_path, monkeypatch)
    expected, expected_log, sample_ids = reference_haplotypes(transcripts, tmp_path, monkeypatch, samples_df,
                                                              colnames, is_x)
    columns, released = fast(transcripts, exons, samples, sample_ids, text, tmp_path)
    cohort = haplotypes.Cohort.from_metadata(samples_df, sample_ids, is_x, PAR1_TO, PAR2_FROM)
    rows, log = fast_haplotypes(transcripts, released, cohort, len(columns))
    assert len(rows) == len(expected)
    for got, (_, want) in zip(rows, expected.iterrows()):
        want = want.tolist()
        if got[1] == 'REF':
            want[5] = want[7] = want[8] = None
        assert got == want
    assert log == expected_log
    # the case reaches what it is meant to: conflicts (removed samples) and merged copies
    assert any(line.split('\t')[1] != '0' for line in log.splitlines()[1:])
    assert expected['Count'].max() > 1
    if n_samples > 100:     # haplotypes with more than 64 copies (frequency strings counted by numpy)
        assert expected[expected['Changes'] != 'REF']['Count'].max() > 64


def test_hash_collision_is_unsupported(monkeypatch):
    """Copies grouped by hash are checked row by row: with weights that cannot tell rows apart, two copies that carry
    different records are caught."""
    np = haplotypes.np
    rows, members = haplotypes.group_copies([np.array([0]), np.array([1])], 2)
    assert sorted(r.tolist() for r in rows) == [[0], [1]]

    class ZeroRng:
        def integers(self, *args, size, **kwargs):
            return np.zeros(size, dtype=np.uint64)
    monkeypatch.setattr(np.random, 'default_rng', lambda seed: ZeroRng())
    with pytest.raises(Unsupported, match='hash'):
        haplotypes.group_copies([np.array([0]), np.array([1])], 2)


def test_metadata_checks():
    samples_df = pd.DataFrame({'Sample name': ['A', 'B', 'B', 'C'], 'Sex': ['male', 'female', 'female', 'male'],
                               'Population code': ['GBR', 'YRI', 'YRI', None],
                               'Superpopulation code': ['EUR', 'AFR', 'AFR', 'EUR']})
    with pytest.raises(Unsupported, match='more than one row'):
        haplotypes.Cohort.from_metadata(samples_df, ['A', 'B'], False, 1, 2)
    with pytest.raises(Unsupported, match='not text'):
        haplotypes.Cohort.from_metadata(samples_df, ['A', 'C'], False, 1, 2)
    with pytest.raises(Unsupported, match='no sample'):
        haplotypes.Cohort.from_metadata(samples_df, [], False, 1, 2)
    with pytest.raises(Unsupported, match='lacks a column'):
        haplotypes.Cohort.from_metadata(samples_df.drop(columns='Sex'), ['A'], False, 1, 2)
    colon = samples_df.assign(**{'Population code': ['GB:R', 'YRI', 'YRI', 'X']})
    with pytest.raises(Unsupported, match="contains ':' or ';'"):
        haplotypes.Cohort.from_metadata(colon, ['A'], False, 1, 2)
