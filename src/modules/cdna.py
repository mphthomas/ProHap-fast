"""The cDNA sequences of some transcripts, as the reference's read_fasta(filename, True) gives them
(legacy/src/modules/common.py), without keeping the whole file in a dict.

read_fasta reads the file line by line: the first line and every later line starting with '>' begin a record, the
record's sequence is every following line without its last character (the newline, or the last base of a last line
that lacks one), and the record is stored under the ID its header gives, a later record replacing an earlier one with
the same ID. Here the file is read in one piece (in text mode, as read_fasta reads it), records are found by '\\n>',
every header is parsed with read_fasta's own rules (below, verbatim), and sequences are built only for the IDs asked
for. A file of ASCII without carriage returns is handled as bytes (the text would be the same); anything else is
decoded as read_fasta decodes it.
"""
import gzip
import re


def element_id(metadata, truncate_accession=True):
    """The ID read_fasta stores a record under; metadata is the header line as read (with its newline)."""
    tag = ''
    accession = ''
    description = ''

    if "|" in metadata:					# the header is at least partially formated
        metadata_parsed = metadata[1:].split('|')
        if 'generic' in metadata_parsed[0]:
            tag = metadata_parsed[0]
        else:
            tag = 'generic_' + metadata_parsed[0]

        if len(metadata_parsed) == 2:					# accession and description are potentially merged -> separate them
            if " " in metadata_parsed[1]:
                accession = metadata_parsed[1].split(' ')[0]
                description = metadata_parsed[1].split(' ', 1)[1]
            else:
                accession = metadata_parsed[1]				# no description -> keep accession as it is
        elif len(metadata_parsed) == 3:					# descripton and accesson already separated -> keep
            accession = metadata_parsed[1]
            description = metadata_parsed[2]

    else:						                    # the header is not formated
        accession = metadata[1:].split(" ")[0]
        if " " in metadata:
            description = metadata.split(" ", 1)[1]

    if (truncate_accession):
        elementID = accession.split('.')[0]
    else:
        elementID = accession
    return elementID


def read_sequences(filename, wanted):
    """{ID: sequence} for the IDs in `wanted` that the file has, as read_fasta(filename, True)[ID]['sequence']."""
    handle = gzip.open(filename, 'rb') if (filename.endswith('.gz')) else open(filename, 'rb')
    data = handle.read()
    handle.close()
    if data.isascii() and b'\r' not in data:
        return _sequences(data, wanted, b'\n', b'\n>', lambda b: b.decode('ascii'))
    # text mode, as read_fasta: universal newlines and the locale's encoding
    handle = gzip.open(filename, 'rt') if (filename.endswith('.gz')) else open(filename, 'r')
    del data
    text = handle.read()
    handle.close()
    return _sequences(text, wanted, '\n', '\n>', lambda t: t)


def _sequences(text, wanted, newline_char, record_start, decode):
    if len(text) == 0:
        return {}
    starts = [0] + [m.start() + 1 for m in re.finditer(re.escape(record_start), text)]
    starts.append(len(text))
    found = {}
    for k in range(len(starts) - 1):
        start, end = starts[k], starts[k + 1]
        newline = text.find(newline_char, start, end)
        if newline == -1:                       # a header that is the file's last line, without a newline
            metadata, body_start = text[start:end], end
        else:
            metadata, body_start = text[start:newline + 1], newline + 1
        eid = element_id(decode(metadata))
        if eid in wanted:
            found[eid] = (body_start, end)
    sequences = {}
    for eid, (body_start, end) in found.items():
        body = text[body_start:end]
        if body and not body.endswith(newline_char):
            body = body[:-1]                    # read_fasta drops the last character of every line
        sequences[eid] = decode(body.replace(newline_char, newline_char[:0]))
    return sequences
