"""The values the reference sees after pandas has read its per-transcript tables (IMPLEMENTATION_NOTES section 3.3).

The reference writes the VCF lines of each transcript into a table with pandas (read_csv, then to_csv) and reads it
back (read_csv again). pandas' C parser reads a table in chunks of rows (low_memory=True, the default) and types each
column of each chunk on its own: a chunk of only numbers becomes int64 or float64, a chunk with any other word stays
text. The chunks are then joined with pandas.core.dtypes.concat.concat_compat. How many rows a chunk has depends on
the width of the table, i.e. on the number of samples (pandas_chunk_rows). So an ID such as `1e5` is `1e5` in one
transcript and `100000.0` in another, depending on the IDs next to it and on the size of the cohort.

emulate_id_column(tokens, chunk_rows) returns, for the ID column of one transcript's table, the strings the reference
gets from str() of each value. It types the chunks itself where the rule is plain (see _type_chunk) and asks pandas
for anything else, and it uses pandas' own concat_compat and to_csv, so it follows pandas rather than a description
of it. tests/test_pandas_emulation.py compares it with a real read_csv -> to_csv -> read_csv round trip.
"""
import re
from functools import lru_cache
from io import StringIO

import numpy as np
import pandas as pd
from pandas.core.dtypes.concat import concat_compat

# pandas' default NA words (pandas._libs.parsers.STR_NA_VALUES)
STR_NA = frozenset(['-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN', '#N/A N/A', '#N/A', 'N/A', 'n/a', 'NA', '<NA>', '#NA',
                    'NULL', 'null', 'NaN', '-NaN', 'nan', '-nan', 'None', ''])
_NUMLIKE = re.compile(r'^\s*[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d*)?\s*$|^\s*[+-]?(inf|infinity|nan)\s*$', re.I)

_INT = re.compile(r'[0-9]{1,18}\Z')
_EXP = re.compile(r'([0-9]{1,15})[eE]([+-]?)([0-9]{1,2})\Z')
_DEC = re.compile(r'[0-9]*\.[0-9]*\Z')
_NAN = float('nan')


def pandas_chunk_rows(n_cols):
    """Rows per chunk when pandas' C parser reads a table n_cols wide: the largest power of two strictly below
    2**20 // n_cols, and at least 1 (pandas._libs.parsers.TextReader, buffer_lines)."""
    heuristic = 2 ** 20 // n_cols
    rows = 1
    while rows * 2 < heuristic:
        rows *= 2
    return rows


def is_safe_str(t):
    """True when pandas certainly keeps t as this exact string, whatever else is in its chunk: not an NA word, not a
    boolean (pandas compares those without case), nothing that might parse as a number, no quote, no surrounding
    whitespace."""
    return not (t in STR_NA or t.lower() in ('true', 'false') or _NUMLIKE.match(t) or '"' in t or t != t.strip())


def _kind(t):
    """('na' | 'int' | 'float' | 'bool' | 'str', value), or (None, None) when the plain rules may not hold."""
    if t in STR_NA:
        return 'na', _NAN
    if is_safe_str(t):
        return 'str', t
    if _INT.match(t):
        return 'int', int(t)
    low = t.lower()
    if low in ('true', 'false'):
        return 'bool', low == 'true'
    m = _EXP.match(t)
    if m:
        exponent = int(m.group(3))
        if exponent <= 22:     # m * 10**e or m / 10**e with both exact doubles: one rounding, as in pandas
            mantissa = float(int(m.group(1)))
            return 'float', (mantissa / 10.0 ** exponent) if m.group(2) == '-' else (mantissa * 10.0 ** exponent)
        return None, None
    if _DEC.match(t) and t != '.' and len(t) <= 16:     # at most 15 digits
        return 'float', float(t)
    return None, None


@lru_cache(maxsize=65536)
def _pandas_chunk(tokens):
    """One chunk typed by pandas itself (two columns, so pandas reads it in one chunk; the second stops a lone empty
    field being taken for a blank line)."""
    text = 'ID\tX\n' + ''.join(t + '\tx\n' for t in tokens)
    return pd.read_csv(StringIO(text), sep='\t')['ID'].to_numpy()


def _type_chunk(tokens):
    """The array pandas' C parser makes of one chunk of a column."""
    kinds = []
    values = []
    for t in tokens:
        k, v = _kind(t)
        if k is None:
            return _pandas_chunk(tuple(tokens))
        kinds.append(k)
        values.append(v)
    present = set(kinds)
    if 'str' in present:
        # text column: every word as it is, NA words missing
        return np.array([_NAN if k == 'na' else t for k, t in zip(kinds, tokens)], dtype=object)
    if 'bool' in present:
        if present == {'bool'}:
            return np.array(values, dtype=bool)
        if present == {'bool', 'na'}:
            return np.array(values, dtype=object)
        return _pandas_chunk(tuple(tokens))
    if present == {'int'}:
        return np.array(values, dtype=np.int64)
    # numbers with NA, or floats: float64. A long integer would be parsed by pandas' float parser, not converted
    if any(k == 'int' and len(t) > 15 for k, t in zip(kinds, tokens)):
        return _pandas_chunk(tuple(tokens))
    return np.array(values, dtype=np.float64)


def read_column(tokens, chunk_rows):
    """The array pandas' C parser makes of a column of a table it reads in chunks of chunk_rows rows."""
    arrays = [_type_chunk(tokens[i:i + chunk_rows]) for i in range(0, len(tokens), chunk_rows)]
    return concat_compat(arrays)


def write_column(values):
    """The text DataFrame.to_csv writes for each value of a column."""
    text = pd.DataFrame({'ID': values, 'X': 'x'}).to_csv(sep='\t', header=False, index=False)
    lines = text.split('\n')
    assert lines[-1] == '' and len(lines) == len(values) + 1
    return [line[:-2] for line in lines[:-1]]


def emulate_id_column(tokens, chunk_rows):
    """str() of each value of the ID column the reference reads back: its VCF tokens through read_csv, to_csv and
    read_csv, all in chunks of chunk_rows rows. Tokens must not contain tabs, newlines or quote characters."""
    if all(is_safe_str(t) for t in tokens):
        return list(tokens)
    first = read_column(list(tokens), chunk_rows)
    second = read_column(write_column(first), chunk_rows)
    return [str(v) for v in second.tolist()]
