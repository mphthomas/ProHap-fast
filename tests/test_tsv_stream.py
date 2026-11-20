"""src/modules/tsv_stream.py, and src/merge_tables.py (D1) against legacy/src/merge_tables.py, as scripts, on crafted
haplotype tables:

    pytest tests/test_tsv_stream.py

pandas reads a table in chunks of rows whose number depends on the table's width (32,768 for the 21 columns of a
haplotype table, 512 for 1,030 columns) and types each column of each chunk on its own. The cases put values whose
text depends on their chunk's type on both sides of chunk edges, at both widths, and check that the fast script writes
the same bytes as the reference, or hand the run to it where joining the chunks rewrites values.
"""
import bz2
import concurrent.futures
import gzip
import io
import itertools
import os
import random
import subprocess
import sys
import warnings
from pathlib import Path

import pytest

pd = pytest.importorskip('pandas')
np = pytest.importorskip('numpy')

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))
from modules.fallback import Unsupported  # noqa: E402
from modules.pandas_emulation import pandas_chunk_rows  # noqa: E402
from modules.tsv_stream import combined_dtype  # noqa: E402

TREES = {'legacy': ROOT / 'legacy' / 'src', 'fast': ROOT / 'src'}
FALLBACK = 'ProHap-fast: falling back to the reference code'
HAPLO_COLUMNS = ['TranscriptID', 'chromosome', 'transcript_biotype', 'HaplotypeID', 'VCF_IDs', 'DNA_changes',
                 'allele_frequencies', 'cDNA_changes', 'all_protein_changes', 'variant_types', 'protein_changes',
                 'reading_frame', 'protein_prefix_length', 'start_missing', 'start_lost', 'splice_sites_affected',
                 'occurrence_count', 'frequency', 'frequency_population', 'frequency_superpopulation', 'samples']
WIDE = 1030             # pandas reads a table this wide in chunks of 512 rows
assert pandas_chunk_rows(len(HAPLO_COLUMNS)) == 32768 and pandas_chunk_rows(WIDE) == 512


# ---------------------------------------------------------------------------------------------------------------
# combined_dtype against pandas


def _array(kind):
    return {
        'int': lambda: np.array([5, 7], dtype=np.int64),
        'float': lambda: np.array([0.5, 1.0]),
        'nan': lambda: np.array([np.nan, np.nan]),
        'bool': lambda: np.array([True, False]),
        'object': lambda: np.array(['a', 3], dtype=object),
        'uint': lambda: np.array([2 ** 63 + 1, 3], dtype=np.uint64),
    }[kind]()


def _text(values):
    return pd.DataFrame({'c': values, 'x': 'x'}).to_csv(sep='\t', header=False, index=False).split('\n')[:-1]


def _frame(kind):
    if kind == 'empty':     # a table with a header and no rows, as read_csv returns it
        return pd.read_csv(io.StringIO('c\tx\n'), sep='\t')
    return pd.DataFrame({'c': _array(kind), 'x': 'x'})


@pytest.mark.parametrize('n_parts', [2, 3])
def test_combined_dtype_follows_pandas(n_parts):
    """For every mix of part dtypes: when combined_dtype accepts it, pandas' concat_compat (the chunks of one file)
    and pd.concat (the tables of several files) write every value as its own part does, and when combined_dtype
    rejects it, at least one of them does not."""
    from pandas.core.dtypes.concat import concat_compat
    kinds = ['int', 'float', 'nan', 'bool', 'object', 'uint', 'empty']
    for combo in itertools.product(kinds, repeat=n_parts):
        parts = [_array(k) for k in combo if k != 'empty']
        own = sum((_text(a) for a in parts), [])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            same_concat = pd.concat([_frame(k) for k in combo])['c'].to_numpy()
            same = _text(same_concat) == own and (not parts or _text(concat_compat(parts)) == own)
        try:
            combined_dtype({a.dtype for a in parts}, 'c')
            accepted = True
        except Unsupported:
            accepted = False
        if accepted:
            assert same, combo
        elif 'empty' not in combo:
            assert not same, combo


def test_float_text_as_object():
    """An object column holds a float chunk's values as Python floats; str() of each must be numpy's text."""
    rng = np.random.default_rng(0)
    values = np.concatenate([rng.random(20000), rng.random(20000) * 10.0 ** rng.integers(-320, 308, 20000),
                             np.frombuffer(rng.bytes(8 * 20000), dtype=np.float64)])
    values = values[np.isfinite(values)]
    assert [str(v) for v in values.tolist()] == values.astype(str).tolist()


# ---------------------------------------------------------------------------------------------------------------
# the scripts


def samples_field(rng, n_samples=40):
    return ';'.join('S%05d:%d' % (rng.randrange(n_samples), rng.randint(1, 2)) for _ in range(rng.randint(1, 6)))


def haplo_rows(rng, n, chromosome, start=0):
    """Rows like those of prohap.py: IDs, changes, frequencies with up to 17 digits, booleans, sample lists."""
    rows = []
    for i in range(start, start + n):
        k = rng.randint(1, 3)
        ids = ';'.join(rng.choice(['rs%d' % rng.randrange(10 ** 6), '%x' % rng.randrange(1, 4096), '2e5', '0012'])
                       for _ in range(k))
        afs = ';'.join(str(round(rng.random() / 3, rng.choice([2, 4, 6]))) for _ in range(k))
        count = rng.randint(1, 500)
        rows.append([
            'ENST%011d' % (i // 7), chromosome, 'protein_coding', 'haplo_chr%s_%x' % (chromosome, i), ids,
            ';'.join('%d:A>G' % rng.randrange(10 ** 7) for _ in range(k)), afs,
            ';'.join('%d:A>G' % rng.randrange(5000) for _ in range(k)), '12:V>12:D', 'SAV', rng.choice(['12:V>12:D', '']),
            str(rng.choice([0, 0, 0, -1])), str(rng.randrange(40)), rng.choice(['False', 'False', 'True']),
            rng.choice(['False', 'True']), rng.choice(['-', '-', '3', '1;2']), str(count), repr(count / 4999),
            'GBR:%r;YRI:0.25' % (rng.random() / 7), 'EUR:%r' % (rng.random() / 3), samples_field(rng)])
    return rows


def text_of(names, rows, end='\n'):
    return '\t'.join(names) + end + ''.join('\t'.join(r) + end for r in rows)


def wide_names():
    return ['TranscriptID', 'chromosome', 'kind', 'HaplotypeID'] + ['c%d' % i for i in range(4, WIDE - 1)] + \
        ['samples']


def wide_rows(rng, n, column_values, start=0):
    """Rows of a WIDE-column table; column_values maps a column index to a function of the row number."""
    rows = []
    for i in range(start, start + n):
        row = ['ENST%d' % (i // 9), '22', 'pc', 'haplo_%x' % i] + [str((i * 7 + c) % 13) for c in range(4, WIDE - 1)]
        row.append(samples_field(rng))
        for c, f in column_values.items():
            row[c] = f(i)
        rows.append(row)
    return rows


def haplo_tables(seed=1, sizes=(300, 450, 0, 120), chromosomes=('1', '21', '22', 'X')):
    rng = random.Random(seed)
    return [text_of(HAPLO_COLUMNS, haplo_rows(rng, n, c)) for n, c in zip(sizes, chromosomes)]


def edge_tables():
    """Two WIDE tables (chunks of 512 rows) whose columns change type at the chunk edges: an object chunk after a
    number or boolean chunk keeps each value's text (`0012` is `12` in the first chunk only), an all-missing chunk
    next to text, a float chunk next to text, and a chunk whose samples lack `:` (typed by pandas)."""
    rng = random.Random(2)
    first = {
        4: lambda i: '0012' if i < 512 or i % 2 else 'abc',      # int chunk, then object chunks
        5: lambda i: 'True' if i < 512 else ['False', 'x'][i % 2],   # bool, then object
        6: lambda i: 'NA' if i < 512 else '1.%d' % i if i < 1024 else 'w',    # all missing, float, object
        7: lambda i: '%r' % (i / 7.0) if i != 1023 else 'z',     # float chunk, object chunk (edge at its last row)
        8: lambda i: 'nan' if i == 511 else 'x' if i == 1050 else str(i),     # float, int, then object chunks
        9: lambda i: 'q' if i == 1050 else str(i),
        WIDE - 1: lambda i: '7' if i == 600 else samples_field(rng),     # the second chunk's samples go through pandas
    }
    one = text_of(wide_names(), wide_rows(rng, 1100, first))
    second = dict(first)
    second[8] = lambda i: 'nan' if i < 512 else 'x'
    second[9] = lambda i: str(i)            # int64 in every chunk of this file; object in the other one
    two = text_of(wide_names(), wide_rows(rng, 700, second))
    return [one, two]


def real_width_edges():
    """A haplotype table of 32,768 + 200 rows (two chunks at its real width), with VCF IDs that are numbers in the
    first chunk and text in the second."""
    rng = random.Random(3)
    rows = haplo_rows(rng, 32968, '22')
    for i, row in enumerate(rows):
        row[4] = '0012' if i < 32768 or i % 3 else 'rs%d' % i
        row[6] = '%r' % (rng.random() / 3)       # every allele frequency a float of up to 17 digits
    return [text_of(HAPLO_COLUMNS, rows)]


def mixed(column, before, after, rows=600):
    """A WIDE table whose column is `before(i)` in the first chunk and `after(i)` in the next."""
    rng = random.Random(4)
    return [text_of(wide_names(), wide_rows(rng, rows, {column: lambda i: before(i) if i < 512 else after(i)}))]


def wide_files(*columns):
    """WIDE tables, one per function, whose column 6 is that function of the row number."""
    rng = random.Random(5)
    return [text_of(wide_names(), wide_rows(rng, 40, {6: f})) for f in columns]


def without_last_newline(tables, table):
    tables[table] = tables[table][:-1]
    return tables


def edit_line(tables, table, line, edit):
    lines = tables[table].split('\n')
    lines[line] = edit(lines[line])
    tables[table] = '\n'.join(lines)
    return tables


def replace_field(tables, table, line, column, value):
    return edit_line(tables, table, line, lambda text: '\t'.join(
        value if i == column else f for i, f in enumerate(text.split('\t'))))


# name: (tables, input suffix, output suffix, fast path?)
CASES = {
    'haplotype_tables': (haplo_tables(), '.tsv.gz', '.tsv.gz', True),
    'plain_files': (haplo_tables(seed=5), '.tsv', '.tsv', True),
    'one_table': (haplo_tables(seed=6, sizes=(500,), chromosomes=('X',)), '.tsv.gz', '.tsv', True),
    'only_headers': (haplo_tables(sizes=(0, 0)), '.tsv.gz', '.tsv.gz', True),
    'no_last_newline': (without_last_newline(haplo_tables(seed=7), 3), '.tsv.gz', '.tsv.gz', True),
    'chunk_edges': (edge_tables(), '.tsv.gz', '.tsv.gz', True),
    'real_width_chunk_edge': (real_width_edges(), '.tsv.gz', '.tsv.gz', True),
    # handed to the reference: pandas joins the chunks by rewriting values
    'int_then_float': (mixed(6, str, lambda i: '%d.5' % i), '.tsv.gz', '.tsv.gz', False),
    'int_then_missing': (mixed(6, str, lambda i: 'NA' if i % 2 else str(i)), '.tsv.gz', '.tsv.gz', False),
    'bool_then_int': (mixed(6, lambda i: 'True', str), '.tsv.gz', '.tsv.gz', False),
    'uint_then_int': (mixed(6, lambda i: str(2 ** 63 + i), str), '.tsv.gz', '.tsv.gz', False),
    'int_file_and_float_file': (wide_files(str, lambda i: '%d.5' % i), '.tsv.gz', '.tsv.gz', False),
    # handed to the reference: input the script does not handle
    'quote': (replace_field(haplo_tables(), 0, 5, 4, '"rs1"'), '.tsv.gz', '.tsv.gz', False),
    'carriage_return': ([t.replace('\n', '\r\n') for t in haplo_tables()], '.tsv.gz', '.tsv.gz', False),
    'non_ascii': (replace_field(haplo_tables(), 1, 3, 2, 'protéin'), '.tsv.gz', '.tsv.gz', False),
    'short_line': (edit_line(haplo_tables(), 1, 7, lambda t: t.rsplit('\t', 1)[0]), '.tsv.gz', '.tsv.gz', False),
    'blank_line': (edit_line(haplo_tables(), 0, 9, lambda t: ''), '.tsv.gz', '.tsv.gz', False),
    'different_columns': ([t.replace('frequency_population', 'pop', 1) if i == 1 else t
                           for i, t in enumerate(haplo_tables())], '.tsv.gz', '.tsv.gz', False),
    'samples_not_last': ([t.replace('\tsamples\n', '\tsamples2\n', 1) for t in haplo_tables()], '.tsv.gz', '.tsv.gz',
                         False),
    'bz2_output': (haplo_tables(), '.tsv.gz', '.tsv.bz2', False),
}


def write(path, text):
    data = text.encode('utf-8')
    if str(path).endswith('.gz'):
        data = gzip.compress(data)
    elif str(path).endswith('.bz2'):
        data = bz2.compress(data)
    path.write_bytes(data)


def read(path):
    data = path.read_bytes()
    if data[:2] == b'\x1f\x8b':
        return gzip.decompress(data)
    if data[:3] == b'BZh':
        return bz2.decompress(data)
    return data


def run_case(case, name):
    tables, suffix, out_suffix, _ = CASES[name]
    inputs = []
    for i, text in enumerate(tables):
        path = case / ('haplo_chr%d%s' % (i, suffix))
        write(path, text)
        inputs.append(str(path))
    result = {}
    env = dict(os.environ, PYTHONHASHSEED='0', PYTHONNOUSERSITE='1')
    for arm, src in TREES.items():
        out = case / arm
        out.mkdir()
        merged = out / ('all' + out_suffix)
        merge = subprocess.run([sys.executable, str(src / 'merge_tables.py'), '-i', ','.join(inputs), '-o',
                                str(merged)], capture_output=True, text=True, env=env)
        result[arm] = merge, (read(merged) if merge.returncode == 0 else None)
    return result


@pytest.fixture(scope='module')
def results(tmp_path_factory):
    """Every case, run by both trees, eight cases at a time."""
    root = tmp_path_factory.mktemp('tables')

    def one(name):
        case = root / name
        case.mkdir()
        return name, run_case(case, name)
    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        return dict(pool.map(one, CASES))


@pytest.mark.parametrize('name', list(CASES))
def test_case(results, name):
    fast_path = CASES[name][3]
    legacy, legacy_output = results[name]['legacy']
    fast, output = results[name]['fast']
    assert legacy.returncode == 0, legacy.stderr
    assert fast.returncode == 0, fast.stderr
    assert (FALLBACK in fast.stderr) == (not fast_path), fast.stderr
    assert output == legacy_output
