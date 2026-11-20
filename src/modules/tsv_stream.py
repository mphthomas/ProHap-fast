"""Stream a tab-separated table through pandas' typing one chunk of rows at a time (IMPLEMENTATION_NOTES section 4).

merge_tables.py and haplo_extract_sample_names.py read whole haplotype tables with pd.read_csv and write them back
with to_csv. That is not a plain copy: pandas' C parser reads a table in chunks of pandas_chunk_rows(width) rows
(low_memory=True, the default) and types each column of each chunk on its own, so a value is written back as pandas
formats its type (`0.17837519623233908` becomes `0.178375196232339`, `True` stays `True`, `0012` becomes `12` in a
chunk of numbers). The fast scripts read the same chunks, run each through read_csv and to_csv on its own, and write
the result straight out, so memory holds one chunk instead of the table.

Typing a chunk on its own gives the reference's text as long as joining the chunks does not change a type. pandas
joins the chunks of one file with concat_compat and the tables of several files with pd.concat: when every chunk of
a column has the same dtype, or one of them is object, every value keeps its text (the object column holds the same
Python values, and str() of a float is numpy's text for it). Any other mix (int64 and float64, bool and a number,
uint64 and int64) rewrites values (`5` becomes `5.0`, `True` becomes `1`), so the script hands the run to the
reference (combined_dtype raises Unsupported). tests/test_tsv_stream.py checks these rules against pandas.

The last column, `samples`, is copied as it is when every value of the chunk contains `:` (`SAMPLE:1;SAMPLE:2`): pandas
keeps such a word as a string, so it is an object column written as it was read. Otherwise that chunk goes through
pandas whole.

Input the scripts do not handle goes to the reference: anything but printable ASCII, tabs and newlines (quote
characters, carriage returns, non-ASCII bytes), a line whose number of fields differs from the header's (pandas skips
blank lines and pads short ones), a header with duplicate or empty names, compression other than gzip.
"""
import io
import re
import subprocess
import warnings
from itertools import chain, islice, repeat

import numpy as np
import pandas as pd
from pandas.io.common import infer_compression

from modules.fallback import Unsupported
from modules.gzip_writer import GzipWriter
from modules.pandas_emulation import pandas_chunk_rows

# tab, newline and printable ASCII other than the quote character: what check_bytes deletes before looking for more
_PLAIN_BYTES = bytes(b for b in range(256) if b in (9, 10) or (0x20 <= b <= 0x7e and b != 0x22))
_OBJECT = np.dtype(object)


def check_compression(path):
    compression = infer_compression(path, 'infer')
    if compression not in (None, 'gzip'):
        raise Unsupported('%s: compressed as %s' % (path, compression))
    return compression == 'gzip'


class Input:
    """The lines of a plain or gzipped file, as bytes; gzip -d decompresses in its own process."""

    def __init__(self, path):
        self.path = path
        self.proc = None
        if check_compression(path):
            self.proc = subprocess.Popen(['gzip', '-dc', path], stdout=subprocess.PIPE, bufsize=1 << 20)
            self.file = self.proc.stdout
        else:
            self.file = open(path, 'rb', buffering=1 << 20)

    def close(self):
        self.file.close()
        if self.proc is not None and self.proc.wait() != 0:
            raise Unsupported('gzip could not decompress ' + self.path)

    def kill(self):
        if self.proc is not None and self.proc.poll() is None:
            self.proc.kill()


class Output:
    """A plain or gzipped (level 1, compressed by a thread) output file, written as bytes."""

    def __init__(self, path):
        self.file = GzipWriter(path) if check_compression(path) else open(path, 'wb', buffering=1 << 20)

    def write(self, data):
        self.file.write(data)

    def write_rows(self, prefixes, lasts):
        """Each prefix, a tab, and the last field with its line end."""
        self.file.write(b''.join(chain.from_iterable(zip(prefixes, repeat(b'\t'), lasts))))

    def close(self):
        self.file.close()


def check_bytes(data, where):
    if data.translate(None, _PLAIN_BYTES):
        raise Unsupported(where + ': a quote, carriage return, control or non-ASCII byte')


def read_header(lines, path):
    """The column names of a table, from its first line."""
    line = next(lines, b'')
    if not line.endswith(b'\n'):
        raise Unsupported(path + ': no header line')
    check_bytes(line, path + ' line 1')
    names = line[:-1].decode('ascii').split('\t')
    if len(set(names)) != len(names) or '' in names:
        raise Unsupported(path + ': duplicate or empty column names')
    if len(names) < 3 or names[-1] != 'samples':
        raise Unsupported(path + ': the last of at least three columns must be samples')
    return names


class Table:
    """One input table, read chunk by chunk as pandas' C parser reads it.

        table = Table(path)
        table.read(write)       # write(prefixes, lasts) for each chunk: the re-formatted text of every column but
                                # the last, and the last column with its line end
        table.dtypes            # per column: the dtypes of its chunks

    A callback rather than a generator, so that no chunk is still referenced while the next one is read.
    """

    def __init__(self, path):
        self.path = path
        self.input = Input(path)
        self.names = read_header(iter(self.input.file), path)
        self.width = len(self.names)
        self.chunk_rows = pandas_chunk_rows(self.width)
        self.dtypes = [set() for _ in self.names]
        self.rows = 0
        self._prefix_header = ('\t'.join(self.names[:-1]) + '\n').encode('ascii')
        self._header = ('\t'.join(self.names) + '\n').encode('ascii')

    def read(self, write):
        lines = iter(self.input.file)
        while True:
            block = list(islice(lines, self.chunk_rows))
            if not block:
                break
            self._chunk(block, write)
            del block
        self.input.close()

    def _chunk(self, block, write):
        """Re-formats one chunk and hands it to write. The last field of each line stays in the line read (a
        memoryview from its last tab), so that the bulk of the text is held once."""
        tabs = self.width - 2       # in a line without its last field
        prefixes = []
        lasts = []
        verbatim = True
        for line in block:
            self.rows += 1
            check_bytes(line, '%s line %d' % (self.path, self.rows + 1))
            k = line.rfind(b'\t')
            prefix = line[:k]
            if k < 0 or prefix.count(b'\t') != tabs:
                raise Unsupported('%s line %d: not %d fields' % (self.path, self.rows + 1, self.width))
            if not line.endswith(b'\n'):
                line += b'\n'
            last = memoryview(line)[k + 1:]
            if verbatim and line.find(b':', k) < 0:
                verbatim = False
            prefixes.append(prefix)
            lasts.append(last)
        if verbatim:
            prefixes = self._reformat(self._prefix_header, prefixes, self.width - 1)
            self.dtypes[-1].add(_OBJECT)
        else:
            # the last column goes through pandas as well
            joined = [p + b'\t' + s[:-1] for p, s in zip(prefixes, lasts)]
            prefixes = lasts = None
            pairs = [line.rpartition(b'\t') for line in self._reformat(self._header, joined, self.width)]
            joined = None
            prefixes = [p for p, _, _ in pairs]
            lasts = [s + b'\n' for _, _, s in pairs]
        block.clear()
        write(prefixes, lasts)

    def _reformat(self, header, lines, n_columns):
        """The lines as pandas reads them as one chunk (their number is at most the chunk size of a table this wide)
        and writes them back; records the dtype of each column."""
        text = header + b'\n'.join(lines) + b'\n'
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            df = pd.read_csv(io.BytesIO(text), sep='\t')
        if df.shape != (len(lines), n_columns):
            raise Unsupported('%s: pandas read a chunk as %d x %d' % (self.path, df.shape[0], df.shape[1]))
        for i, dtype in enumerate(df.dtypes):
            self.dtypes[i].add(dtype)
        out = df.to_csv(sep='\t', header=False, index=False).encode('ascii').split(b'\n')
        if out[-1] != b'' or len(out) != len(lines) + 1:
            raise Unsupported(self.path + ': pandas wrote a chunk on other than one line per row')
        out.pop()
        return out

    def kill(self):
        self.input.kill()


def combined_dtype(dtypes, column):
    """The dtype pandas gives a column joined from parts of these dtypes (None for no parts), when every value keeps
    its text; raises Unsupported when joining them would rewrite values."""
    if not dtypes:
        return None
    if len(dtypes) == 1:
        return next(iter(dtypes))
    if _OBJECT in dtypes:
        return _OBJECT
    raise Unsupported('column %s has chunks typed %s, which pandas joins by rewriting values'
                      % (column, ' and '.join(sorted(str(d) for d in dtypes))))
