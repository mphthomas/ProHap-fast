"""src/modules/pandas_emulation.py against pandas itself:

    pytest tests/test_pandas_emulation.py

The reference reads each transcript's VCF lines through pandas (read_csv, to_csv, read_csv); pandas types each chunk of
rows on its own, and the chunk size depends on the number of columns. The round-trip tests build tables as wide as a
cohort of about 1,000, 16,000, 65,000 and 350,000 samples, so that pandas reads them in chunks of 512, 32, 8 and 1
rows, and compare every ID with emulate_id_column.
"""
import io
import math
import random
import sys
import warnings
from pathlib import Path

import pytest

pd = pytest.importorskip('pandas')
np = pytest.importorskip('numpy')

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))
from modules import pandas_emulation as pe  # noqa: E402

# what the IDs of a VCF can look like, weighted towards tokens pandas retypes: hex line numbers (what `.` becomes),
# numbers, NA words, booleans, rs IDs and words
FIXED = ['.', 'rs1', 'rs123', 'abc', 'e5', '1e', '1_000', '0x1A', 'NA', 'nan', 'NaN', 'null', 'None', 'N/A', 'True',
         'False', 'TRUE', 'true', 'tRuE', 'FALSE', '0', '00', '0012', '7', '123', '1.50', '1.', '.5', '2e5', '1E5',
         '1e+5', '1e-5', '3e22', '3e23', '1e400', 'inf', '-inf', 'Infinity', '+5', '-1', '-0', '1d5',
         '123456789012345', '1234567890123456', '123456789012345678', '1234567890123456789',
         '12345678901234567890', '0.1', '0.30000000000000004', '1.5e3', 'X1', '1;2', 'rs1;rs2', 'ID:1']


def random_token(rng):
    r = rng.random()
    if r < 0.35:
        return rng.choice(FIXED)
    if r < 0.75:
        return '%x' % rng.randint(1, 1 << rng.choice([8, 12, 16, 20, 24]))      # hex line numbers
    if r < 0.85:
        return str(rng.randint(0, 10 ** rng.randint(1, 15)))
    if r < 0.92:
        return '%d.%d' % (rng.randint(0, 9999), rng.randint(0, 99999))
    return 'rs%d' % rng.randint(1, 10 ** 9)


def same_value(a, b):
    if type(a) is not type(b):
        return False
    if isinstance(a, float) and math.isnan(a):
        return math.isnan(b)
    return a == b


def pandas_chunk(tokens):
    text = 'ID\tX\n' + ''.join(t + '\tx\n' for t in tokens)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return pd.read_csv(io.StringIO(text), sep='\t')['ID'].to_numpy()


# ---------------------------------------------------------------------------------------------------------------
# pieces


@pytest.mark.parametrize('n_cols, rows', [(9 + 1016, 512), (2048, 256), (16384, 32), (31775, 32), (31776, 16),
                                          (65536, 8), (116508, 8), (116509, 4), (349525, 2), (349526, 1)])
def test_pandas_chunk_rows(n_cols, rows):
    assert pe.pandas_chunk_rows(n_cols) == rows


@pytest.mark.parametrize('n_cols', [1025, 65536])
def test_chunk_rows_is_what_pandas_does(n_cols):
    """A word at the start of the second chunk turns that chunk, and only that chunk, into text."""
    b = pe.pandas_chunk_rows(n_cols)
    ids = ['5'] * (3 * b)
    ids[b] = 'abc'
    filler = '\t'.join('x' for _ in range(n_cols - 1))
    text = '\t'.join(['ID'] + ['c%d' % i for i in range(n_cols - 1)]) + '\n' + ''.join(t + '\t' + filler + '\n'
                                                                                         for t in ids)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        column = pd.read_csv(io.StringIO(text), sep='\t')['ID'].tolist()
    assert [isinstance(v, str) for v in column] == [False] * b + [True] * b + [False] * b


def test_is_safe_str():
    for t in ['rs1', '.', 'abc', 'e5', '1_000', '0x1A', 'X1', 'ID:1', '1;2']:
        assert pe.is_safe_str(t), t
    for t in ['NA', 'nan', '', 'True', 'tRuE', '0', '1e', '0012', '2e5', '1.', '.5', 'inf', '-Infinity', '+5', '1d5',
              ' x', 'x ', 'a"b']:
        assert not pe.is_safe_str(t), t


def test_float_tokens_parse_as_in_pandas():
    """The words _kind turns into floats itself get the same double as pandas' parser gives them."""
    rng = random.Random(1)
    tokens = []
    for _ in range(20000):
        digits = rng.randint(1, 15)
        mantissa = str(rng.randint(0, 10 ** digits - 1)).zfill(rng.choice([1, digits]))
        r = rng.random()
        if r < 0.4:
            tokens.append(mantissa + rng.choice('eE') + rng.choice(['', '+', '-']) + str(rng.randint(0, 22)))
        else:
            cut = rng.randint(0, len(mantissa))
            tokens.append(mantissa[:cut] + '.' + mantissa[cut:])
    tokens = [t for t in tokens if pe._kind(t)[0] == 'float']
    assert len(tokens) > 19000
    parsed = pandas_chunk(tokens)
    assert parsed.dtype == np.float64
    ours = np.array([pe._kind(t)[1] for t in tokens])
    assert (parsed.view(np.int64) == ours.view(np.int64)).all()


def test_type_chunk_matches_pandas():
    rng = random.Random(2)
    for trial in range(3000):
        n = rng.choice([1, 1, 2, 3, 5, 8, 16, 32])
        pool = rng.sample(FIXED, 6) + ['%x' % rng.randint(1, 4096) for _ in range(4)]
        tokens = [rng.choice(pool) for _ in range(n)]
        expected = pandas_chunk(tokens)
        got = pe._type_chunk(tokens)
        assert got.dtype == expected.dtype, tokens
        assert all(same_value(a, b) for a, b in zip(got.tolist(), expected.tolist())), tokens


# ---------------------------------------------------------------------------------------------------------------
# the whole round trip, as the reference does it


def reference_ids(ids, n_cols, tmp_path, usecols=False):
    """str() of each ID as the reference sees it: parse_vcf's read_csv -> to_csv, get_haplotypes' read_csv. With
    usecols, pandas converts only POS and ID (much faster at 350,000 columns; it still tokenizes, and chunks, the
    whole width: test_usecols_keeps_the_chunks), and the file read back is the to_csv text of those two columns with
    the other columns' text appended."""
    filler = '\t'.join('0|0' for _ in range(n_cols - 2))
    header = '\t'.join(['POS', 'ID'] + ['S%d' % i for i in range(n_cols - 2)]) + '\n'
    lines = ''.join('%d\t%s\t%s\n' % (1000 + i, t, filler) for i, t in enumerate(ids))
    path = tmp_path / 'transcript.tsv'
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        if usecols:
            first = pd.read_csv(io.StringIO(header + lines), sep='\t', usecols=['POS', 'ID'])
            text = first.to_csv(sep='\t', index=False, header=False)
            path.write_text(header + ''.join(line + '\t' + filler + '\n' for line in text.splitlines()))
            column = pd.read_csv(path, sep='\t', usecols=['ID'])['ID'].tolist()
        else:
            pd.read_csv(io.StringIO(header + lines), sep='\t').to_csv(path, sep='\t', index=False, header=True)
            column = pd.read_csv(path, sep='\t')['ID'].tolist()
    return [str(v) for v in column]


def test_usecols_keeps_the_chunks():
    n_cols = 350000
    b = pe.pandas_chunk_rows(n_cols)
    assert b == 1
    filler = '\t'.join('x' for _ in range(n_cols - 1))
    text = '\t'.join(['ID'] + ['c%d' % i for i in range(n_cols - 1)]) + '\n' + ''.join(t + '\t' + filler + '\n'
                                                                                         for t in ['5', 'abc', '7'])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        column = pd.read_csv(io.StringIO(text), sep='\t', usecols=['ID'])['ID'].tolist()
    assert column == [5, 'abc', 7]


@pytest.mark.parametrize('n_cols, n_rows, seeds, usecols', [(1025, 1200, 4, False), (16384, 160, 2, False),
                                                            (65536, 64, 3, True), (350000, 40, 3, True)],
                         ids=['chunk512', 'chunk32', 'chunk8', 'chunk1'])
def test_emulate_id_column_round_trip(n_cols, n_rows, seeds, usecols, tmp_path):
    b = pe.pandas_chunk_rows(n_cols)
    for seed in range(seeds):
        rng = random.Random(seed)
        # some tables mostly of numbers, so that whole chunks are numeric
        numeric = seed % 2 == 0
        ids = []
        for _ in range(n_rows):
            if numeric and rng.random() < 0.8:
                ids.append(rng.choice(['%x' % rng.randint(1, 255), str(rng.randint(0, 999)), '2e5', 'NA', 'nan',
                                       '1.50', 'True', 'False', '0012', '1e3']))
            else:
                ids.append(random_token(rng))
        expected = reference_ids(ids, n_cols, tmp_path, usecols)
        assert pe.emulate_id_column(ids, b) == expected, (n_cols, seed)
        if numeric:
            assert expected != ids          # pandas did retype something


def test_all_safe_tokens_are_kept():
    ids = ['rs1', 'abc', '.', 'X1']
    assert pe.emulate_id_column(ids, 1) == ids
