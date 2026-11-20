"""Reading a FASTA file record by record, exactly as modules.common.read_fasta reads it (IMPLEMENTATION_NOTES 8.3).

read_fasta returns a dict keyed by accession: a record whose accession was seen before replaces the earlier one's
value but keeps its position. `records` yields the records in file order, with the same parsing, and `Accessions`
tells afterwards whether any accession repeated (then a caller that streamed must redo the work the dict way).
"""
import gzip
from array import array

import numpy as np


def records(filename, truncate_accession=False):
    """(element ID, {'tag', 'accession', 'description', 'sequence'}) for every record, in file order: read_fasta's
    loop, yielding instead of storing."""
    fasta_file = gzip.open(filename, 'rt') if (filename.endswith('.gz')) else open(filename, 'r')
    try:
        metadata = fasta_file.readline()
        parts = []
        while metadata != "":
            line = fasta_file.readline()
            while not line.startswith('>') and not line == "":
                parts.append(line[:-1])
                line = fasta_file.readline()
            sequence = ''.join(parts)
            parts = []

            tag = ''
            accession = ''
            description = ''

            if "|" in metadata:
                metadata_parsed = metadata[1:].split('|')
                if 'generic' in metadata_parsed[0]:
                    tag = metadata_parsed[0]
                else:
                    tag = 'generic_' + metadata_parsed[0]

                if len(metadata_parsed) == 2:
                    if " " in metadata_parsed[1]:
                        accession = metadata_parsed[1].split(' ')[0]
                        description = metadata_parsed[1].split(' ', 1)[1]
                    else:
                        accession = metadata_parsed[1]
                elif len(metadata_parsed) == 3:
                    accession = metadata_parsed[1]
                    description = metadata_parsed[2]

            else:
                accession = metadata[1:].split(" ")[0]
                if " " in metadata:
                    description = metadata.split(" ", 1)[1]

            element_id = accession.split('.')[0] if truncate_accession else accession
            yield element_id, {'tag': tag, 'accession': accession, 'description': description.replace('\n', ''),
                               'sequence': sequence}
            metadata = line
    finally:
        fasta_file.close()


class Accessions:
    """64-bit hashes of the element IDs seen, to find out at the end whether one repeated. Equal hashes of different
    IDs count as a repeat (the caller then takes the exact, slower way)."""

    def __init__(self):
        self.hashes = array('q')

    def add(self, element_id):
        self.hashes.append(hash(element_id))

    def repeated(self):
        h = np.frombuffer(self.hashes, dtype=np.int64) if len(self.hashes) else np.empty(0, dtype=np.int64)
        return len(np.unique(h)) != len(h)

    def __len__(self):
        return len(self.hashes)
