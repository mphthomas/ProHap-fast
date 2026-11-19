"""Assigning VCF records to transcripts (IMPLEMENTATION_NOTES section 3.2): the reference's parse_vcf without its
per-transcript tables.

The reference (legacy/src/modules/vcf_reader.py) sweeps the VCF once, keeping a queue of the transcripts that may
still receive records, and writes every transcript's lines to a table with one column per sample. Here the same
sweep, step for step (the same queue, the same conditions for reading on and for taking transcripts off the queue,
the same validity check, exon test and line count for `.` IDs), yields each transcript with the records it received,
in the order the reference's queue releases them. A record is parsed once, the first time it falls in an exon, and
shared by every transcript it belongs to; of its genotypes only the haplotype copies that carry the ALT allele are
kept (`copies`: 2 * rank for the first copy of the rank-th included sample, 2 * rank + 1 for the second).

Anything the fast path does not reproduce exactly raises Unsupported, and prohap.py hands the run to the reference:
a carriage return or non-ASCII byte, a header that does not start with the nine VCF columns, an empty or repeated
column name, a quote character, whitespace other than single tabs between the first nine fields, a POS that is not
plain digits, FORMAT other than GT, INFO that pandas might read as something other than text, genotypes other than
`a|b` with a and b one of 0-9 or `.` (unphased `a/b` included), or a record with a different number of samples.
"""
import bisect
import re
import subprocess

import numpy as np

from modules.fallback import Unsupported
from modules.legacy_import import legacy
from modules.pandas_emulation import emulate_id_column, is_safe_str

FIXED_COLUMNS = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
_POS = re.compile(r'(0|[1-9][0-9]{0,17})\Z')
REF_WORD = np.frombuffer(b'0|0\t', dtype='<u4')[0]
_IS_ALLELE = np.zeros(256, dtype=bool)
_IS_ALLELE[[ord(c) for c in '0123456789.']] = True
_TAB, _NL, _PIPE, _ONE = (ord(c) for c in '\t\n|1')


class VCFReader:
    """The lines of a VCF as bytes; a .gz file is decompressed by a gzip process."""

    def __init__(self, path):
        self.process = None
        if path.endswith('.gz'):
            self.process = subprocess.Popen(['gzip', '-dc', path], stdout=subprocess.PIPE, bufsize=1 << 20)
            self.handle = self.process.stdout
        else:
            self.handle = open(path, 'rb', buffering=1 << 20)
        self.path = path

    def readline(self):
        line = self.handle.readline()
        if b'\r' in line or not line.isascii():
            raise Unsupported('carriage return or non-ASCII byte in ' + self.path)
        return line

    def close(self):
        self.handle.close()
        if self.process is not None:
            if self.process.poll() is None:
                self.process.kill()
            self.process.wait()

    def finish(self):
        """Close after reading everything; a gzip error means the reference would have stopped too."""
        self.handle.close()
        if self.process is not None and self.process.wait() != 0:
            raise Unsupported('gzip could not decompress ' + self.path)


def read_header(reader):
    """As parse_vcf: the header is every line starting with '#' before the first other line, and the last of them
    names the columns. Returns (line number of the first data line, column names, first data line); the data line is
    b'' when there is none, and then nothing else is checked (the reference writes empty output)."""
    line_number = 1
    header = None
    line = reader.readline()
    while line != b'' and line.startswith(b'#'):
        header = line[1:]
        line_number += 1
        line = reader.readline()
    if line == b'':
        return line_number, None, line
    if header is None or not header.endswith(b'\n'):
        raise Unsupported('no header line naming the columns')
    columns = header[:-1].decode('ascii').split('\t')
    samples = columns[9:]
    if columns[:9] != FIXED_COLUMNS:
        raise Unsupported('the header does not start with the nine VCF columns')
    if not samples:
        raise Unsupported('the VCF has no samples')
    if any(s == '' or s != s.strip() or '"' in s or not s.isprintable() for s in samples):
        raise Unsupported('an empty sample name, or one with surrounding whitespace or quotes')
    if len(set(columns)) != len(columns):
        raise Unsupported('a column name is repeated in the header')
    return line_number, columns, line


def exon_index(exons):
    """(starts, prefix maximum of ends) of exons ordered by start: a record overlaps an exon in the reference's sense
    (exon.start < pos + len(REF) and exon.end >= pos, scanning exons by start) iff k = bisect_left(starts,
    pos + len(REF)) > 0 and max_ends[k - 1] >= pos."""
    starts = [e[0] for e in exons]
    max_ends = []
    top = None
    for _, end in exons:
        top = end if top is None or end > top else top
        max_ends.append(top)
    return starts, max_ends


class RecordParser:
    """Turns data lines into records: (POS text, POS, ID, REF, ALT, INFO, copies). n_vcf_samples is the number of
    sample columns; rank[c] is the position of sample column c among the included samples, or -1."""

    def __init__(self, n_vcf_samples, rank):
        self.n = n_vcf_samples
        self.rank = np.asarray(rank, dtype=np.int64)
        self.valid = {}

    def head(self, line):
        """(tabs, fields) of the first nine fields of a line, checked."""
        tabs = []
        t = -1
        for _ in range(9):
            t = line.find(b'\t', t + 1)
            if t < 0:
                raise Unsupported('a data line with fewer than ten fields')
            tabs.append(t)
        head = line[:tabs[8]].decode('ascii')
        if len(head.split()) != 9 or head.count('\t') != 8 or '"' in head or not head.replace('\t', ' ').isprintable():
            raise Unsupported('whitespace other than single tabs, a quote or a control character in: ' + head[:200])
        fields = head.split('\t')
        if not _POS.match(fields[1]):
            raise Unsupported('POS is not plain digits: ' + fields[1])
        if fields[8] != 'GT':
            raise Unsupported('FORMAT is not GT: ' + fields[8])
        return tabs, fields

    def is_valid(self, ref, alt):
        """legacy check_vcf_line_validity with min_af 0 (prohap.py hands any other -af to the reference)."""
        key = (ref, alt)
        valid = self.valid.get(key)
        if valid is None:
            valid = self.valid[key] = legacy.vcf_reader.check_vcf_line_validity('', 0, ref, alt)
        return valid

    def copies(self, line, start):
        """The haplotype copies that carry the ALT allele, from the genotypes starting at byte `start`."""
        if not line.endswith(b'\n'):
            line += b'\n'
        if len(line) - start != 4 * self.n:
            raise Unsupported('a record whose genotypes are not %d tokens a|b' % self.n)
        words = np.frombuffer(line, dtype='<u4', offset=start)
        other = np.flatnonzero(words != REF_WORD)
        w = words[other]
        a = w & 0xff
        b = (w >> 16) & 0xff
        end = np.where(other == self.n - 1, _NL, _TAB)
        if not (_IS_ALLELE[a] & (((w >> 8) & 0xff) == _PIPE) & _IS_ALLELE[b] & ((w >> 24) == end)).all():
            raise Unsupported('a genotype other than a|b (one character per allele, phased)')
        first = self.rank[other[a == _ONE]]
        second = self.rank[other[b == _ONE]]
        copies = np.concatenate((2 * first[first >= 0], 2 * second[second >= 0] + 1))
        copies.sort()
        return copies.astype(np.int32)


def assign_records(reader, first_line, first_line_number, transcripts, parser):
    """The reference's sweep. transcripts: objects with start, end, exon_starts and exon_max_ends, in the order of
    all_transcripts (sorted by start). Yields (index of the transcript, its records in VCF order) as the reference's
    queue releases them."""
    head = parser.head(first_line)
    state = {'line': first_line, 'head': head, 'pos': int(head[1][1]), 'count': first_line_number}
    queue = []          # [end, index, records], ordered by end; ties: the newer transcript first (bisect_left)
    queue_ends = []

    def read_until(limit, finalize):
        line, head, pos, count = state['line'], state['head'], state['pos'], state['count']
        while (finalize or pos < limit) and line != b'':
            tabs, fields = head
            ref = fields[3]
            alt = fields[4]
            if parser.is_valid(ref, alt):
                record = None
                end = pos + len(ref)
                for entry in queue:
                    t = transcripts[entry[1]]
                    k = bisect.bisect_left(t.exon_starts, end)
                    if k and t.exon_max_ends[k - 1] >= pos:
                        if record is None:
                            info = fields[7]
                            if not is_safe_str(info):
                                raise Unsupported('INFO that pandas might not read as text: ' + info[:200])
                            vcf_id = fields[2] if fields[2] != '.' else hex(count)[2:]
                            record = (fields[1], pos, vcf_id, ref, alt, info, parser.copies(line, tabs[8] + 1))
                        entry[2].append(record)
            count += 1
            line = reader.readline()
            if line == b'':
                break
            head = parser.head(line)
            pos = int(head[1][1])
        state.update(line=line, head=head, pos=pos, count=count)

    def release(finalize):
        while queue and (finalize or queue[0][0] < state['pos']):
            entry = queue.pop(0)
            queue_ends.pop(0)
            yield entry[1], entry[2]

    last_start = None
    for index, t in enumerate(transcripts):
        if last_start is not None and last_start < t.start:
            read_until(t.start, False)
            yield from release(False)
        i = bisect.bisect_left(queue_ends, t.end)
        queue.insert(i, [t.end, index, []])
        queue_ends.insert(i, t.end)
        last_start = t.start
    read_until(None, True)
    yield from release(True)


def transcript_columns(records, chunk_rows):
    """The columns get_haplotypes reads back from a transcript's table, as str() of each value: POS, ID (through
    pandas_emulation), REF, ALT, INFO. POS, REF, ALT and INFO are kept as text by pandas (checked by RecordParser)."""
    return {
        'POS': [r[0] for r in records],
        'ID': emulate_id_column([r[2] for r in records], chunk_rows),
        'REF': [r[3] for r in records],
        'ALT': [r[4] for r in records],
        'INFO': [r[5] for r in records],
    }
