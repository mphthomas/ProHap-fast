"""src/modules/tsv_stream.py: the rule by which pandas joins chunks and tables without rewriting values.

    pytest tests/test_tsv_stream.py
"""
import io
import itertools
import sys
import warnings
from pathlib import Path

import pytest

pd = pytest.importorskip('pandas')
np = pytest.importorskip('numpy')

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))
from modules.fallback import Unsupported  # noqa: E402
from modules.tsv_stream import combined_dtype  # noqa: E402


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
