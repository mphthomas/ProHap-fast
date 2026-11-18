"""Byte-for-byte comparison of two output files after decompression.

Gzip headers carry a timestamp, so two runs of the same code write different .gz files; the contract is about the
decompressed bytes. A file is read as gzip when it starts with the gzip magic bytes, whatever its name, and a gzip
file never equals a plain one, even with the same content: the next step would read them differently.
"""
import gzip
import hashlib
import itertools
import os

GZIP_MAGIC = b'\x1f\x8b'
IDENTICAL, DIFFERENT, MISSING, ABSENT = 'IDENTICAL', 'DIFFERENT', 'MISSING', 'ABSENT'
SNIPPET = 100


class Result:
    def __init__(self, name, status, sha_a=None, sha_b=None, line=None, text_a=None, text_b=None, note=''):
        self.name = name
        self.status = status      # IDENTICAL, DIFFERENT, MISSING (in one arm or, for a required file, in both)
        self.sha_a = sha_a        # sha256 of the decompressed bytes; None when the file is missing
        self.sha_b = sha_b
        self.line = line          # first differing line, counted from 1
        self.text_a = text_a      # that line in each file, shortened; None at the end of the file
        self.text_b = text_b
        self.note = note

    @property
    def ok(self):
        return self.status in (IDENTICAL, ABSENT)

    def difference(self):
        """Short description of the first difference, '' when there is none."""
        if self.note:
            return self.note
        if self.line is not None:
            return 'line %d' % self.line
        return ''


def is_gzip(path):
    with open(path, 'rb') as fh:
        return fh.read(2) == GZIP_MAGIC


def _open(path):
    return gzip.open(path, 'rb') if is_gzip(path) else open(path, 'rb')


def _show(line):
    if line is None:
        return None
    text = line.rstrip(b'\n').decode('utf-8', 'replace')
    return text if len(text) <= SNIPPET else text[:SNIPPET] + '...'


def sha256(path):
    """sha256 of a file's decompressed bytes."""
    h = hashlib.sha256()
    with _open(path) as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def compare(name, a, b, required=True):
    """Compare file a (legacy) with file b (fast). name labels the result. A file missing from both arms is MISSING
    when required, else ABSENT (a file the workflow writes only in some cases)."""
    have_a, have_b = os.path.isfile(a), os.path.isfile(b)
    if not (have_a or have_b):
        return Result(name, MISSING if required else ABSENT, note='in neither arm' if required else '')
    if not (have_a and have_b):
        present = a if have_a else b
        sha = sha256(present)
        return Result(name, MISSING, sha if have_a else None, None if have_a else sha,
                      note='only in the %s arm' % ('legacy' if have_a else 'fast'))
    try:
        gz_a, gz_b = is_gzip(a), is_gzip(b)
        ha, hb = hashlib.sha256(), hashlib.sha256()
        line = text_a = text_b = None
        with _open(a) as fa, _open(b) as fb:
            for n, (la, lb) in enumerate(itertools.zip_longest(fa, fb), 1):
                if la is not None:
                    ha.update(la)
                if lb is not None:
                    hb.update(lb)
                if line is None and la != lb:
                    line, text_a, text_b = n, _show(la), _show(lb)
    except (OSError, EOFError, gzip.BadGzipFile) as e:
        return Result(name, DIFFERENT, note='cannot read: %s' % e)
    sha_a, sha_b = ha.hexdigest(), hb.hexdigest()
    if gz_a != gz_b:
        return Result(name, DIFFERENT, sha_a, sha_b, line, text_a, text_b,
                      note='%s is gzip, %s is not' % (('legacy', 'fast') if gz_a else ('fast', 'legacy')))
    if line is None:
        return Result(name, IDENTICAL, sha_a, sha_b)
    return Result(name, DIFFERENT, sha_a, sha_b, line, text_a, text_b)
