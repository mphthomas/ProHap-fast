import argparse
import os
import re
from array import array
from itertools import islice

import gffutils
import numpy as np
import pandas as pd

from modules.fallback import Unsupported, to_legacy
from modules.gzip_writer import GzipWriter
from modules.pandas_emulation import is_safe_str
from modules.tsv_stream import Input, Output, check_bytes, check_compression

# ProHap-fast (IMPLEMENTATION_NOTES section 4.4): the same arguments and outputs as the reference
# (legacy/src/fasta_simplify_headers.py), which reads the whole FASTA into a dict, the whole haplotype table into
# pandas, and collects every header row before writing them with to_csv. This script streams the FASTA record by
# record, keeps only HaplotypeID -> TranscriptID from the table (sorted fixed-width byte strings and transcript codes,
# about 25 bytes per row), caches each transcript's gene names, and writes the header table with to_csv one block of
# rows at a time. Input it does not handle goes to the reference: control or non-ASCII bytes or carriage returns
# anywhere, quotes in the table, two records with the same accession (the reference's dict would keep the first
# one's place and the last one's content), a table pandas would read differently from its lines, a looked-up
# HaplotypeID or TranscriptID that pandas might not keep as that text, and every input on which the reference stops.
# Gzipped outputs are written at level 1 (the reference: 9 for the header table, gzip.open's 9 for the FASTA).

parser = argparse.ArgumentParser(description='Reads a FASTA file, extracts all extra information in headers into a separate tab-separated file.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input FASTA file", metavar="FILE")

parser.add_argument("-hap_tsv", dest="haplo_table", required=False,
                    help="input haplotype table", metavar="FILE")

parser.add_argument("-hap_prefix", dest="haplo_prefix", required=False,
                    help="prefix for haplotype protein ID (default: 'haplo_')", default='haplo_')

parser.add_argument("-var_prefix", dest="var_prefix", required=False,
                    help="prefix for variant protein ID (default: 'var_')", default='var_')

parser.add_argument("-db", dest="annotation_db", required=True,
                    help="DB file created by gffutils from GTF")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output FASTA file", metavar="FILE")

parser.add_argument("-header", dest="output_file_header", required=True,
                    help="output header file (tab-separated)", metavar="FILE")

args = parser.parse_args()

HEADER_COLUMNS = ['accession', 'tag', 'matching_proteins', 'position_within_protein', 'start', 'reading_frame']
BLOCK = 1 << 14             # FASTA records per block
TABLE_CHUNK = 1 << 20       # table rows per array of IDs while reading
MAX_ID_BYTES = 4 << 30      # the sorted IDs (rows x longest ID) may take this much
# tab, newline and printable ASCII: what a FASTA line may hold (the reference reads it in text mode, which would
# turn a carriage return into a line end)
_FASTA_BYTES = bytes(b for b in range(256) if b in (9, 10) or 0x20 <= b <= 0x7e)


class HaplotypeTable:
    """HaplotypeID -> TranscriptID of a haplotype table, as `pd.read_table(...).set_index('HaplotypeID')` followed by
    `.loc[ID]['TranscriptID']` gives them.

    Every line is checked as tsv_stream checks a table (pandas must read it as its lines say). Lookups are exact only
    when the index is unique text: every HaplotypeID must be a word that pandas keeps as text whatever else is in its
    chunk (is_safe_str), and no two may be equal; otherwise the first lookup hands the run to the reference (which
    never looks anything up when no record needs it, so such a table is fine until then). A TranscriptID that is
    looked up must be such a word too, so that the value is its text in any chunk."""

    def __init__(self, path):
        self.problem = None
        source = Input(path)
        try:
            self._read(source, path)
        except Unsupported:
            source.kill()
            raise
        source.close()

    def _read(self, source, path):
        lines = iter(source.file)
        header = next(lines, b'')
        if not header.endswith(b'\n'):
            raise Unsupported(path + ': no header line')
        check_bytes(header, path + ' line 1')
        names = header[:-1].decode('ascii').split('\t')
        if len(set(names)) != len(names) or '' in names:
            raise Unsupported(path + ': duplicate or empty column names')
        if 'HaplotypeID' not in names:
            raise Unsupported(path + ': no HaplotypeID column')    # set_index stops the reference
        id_col = names.index('HaplotypeID')
        tr_col = names.index('TranscriptID') if 'TranscriptID' in names else None
        if tr_col is None:
            self.problem = path + ': no TranscriptID column'
        tabs = len(names) - 1
        split = max(id_col, tr_col or 0) + 1
        codes = array('i')
        code_of = {}
        self.transcripts = []
        parts = []
        ids = []
        row = 0
        for line in lines:
            row += 1
            check_bytes(line, '%s line %d' % (path, row + 1))
            if line.count(b'\t') != tabs:
                raise Unsupported('%s line %d: not %d fields' % (path, row + 1, tabs + 1))
            if line.endswith(b'\n'):
                line = line[:-1]
            fields = line.split(b'\t', split)
            haplotype = fields[id_col]
            if self.problem is None and not is_safe_str(haplotype.decode('ascii')):
                self.problem = '%s line %d: HaplotypeID %r, which pandas may read as other than this text' % (
                    path, row + 1, haplotype.decode('ascii'))
            ids.append(haplotype)
            if tr_col is not None:
                transcript = fields[tr_col]
                code = code_of.get(transcript)
                if code is None:
                    code = code_of[transcript] = len(self.transcripts)
                    self.transcripts.append(transcript.decode('ascii'))
                codes.append(code)
            if len(ids) == TABLE_CHUNK:
                parts.append(np.array(ids, dtype=np.bytes_))
                ids = []
        parts.append(np.array(ids, dtype=np.bytes_))
        del ids
        width = max(p.dtype.itemsize for p in parts)
        if width * row > MAX_ID_BYTES:
            raise Unsupported('%s: haplotype IDs of up to %d bytes in %d rows' % (path, width, row))
        keys = np.concatenate(parts).astype('S%d' % width)
        del parts
        order = np.argsort(keys)
        self.keys = keys[order]
        del keys
        self.codes = np.frombuffer(codes, dtype=np.int32)[order] if tr_col is not None else None
        del order, codes
        # Byte strings hold no NUL here (check_bytes), so numpy's comparison is exact.
        if self.problem is None and row > 1 and bool((self.keys[1:] == self.keys[:-1]).any()):
            self.problem = path + ': a HaplotypeID on more than one row'
        self.safe = {}      # transcript code -> its text is kept by pandas

    def transcripts_of(self, haplotypes):
        """The TranscriptID of each haplotype ID, in order."""
        if self.problem is not None:
            raise Unsupported(self.problem)
        width = self.keys.dtype.itemsize
        for haplotype in haplotypes:
            if not is_safe_str(haplotype) or len(haplotype) > width:
                raise Unsupported('haplotype ID %r is not in the table, or pandas may not read it as text'
                                  % haplotype)
        wanted = np.array(haplotypes, dtype=self.keys.dtype)
        at = np.searchsorted(self.keys, wanted)
        found = at < len(self.keys)
        found[found] = self.keys[at[found]] == wanted[found]
        if not found.all():
            raise Unsupported('haplotype ID %r is not in the table' % haplotypes[int(np.argmin(found))])
        result = []
        for code in self.codes[at].tolist():
            safe = self.safe.get(code)
            if safe is None:
                safe = self.safe[code] = is_safe_str(self.transcripts[code])
            if not safe:
                raise Unsupported('TranscriptID %r, which pandas may read as other than this text'
                                  % self.transcripts[code])
            result.append(self.transcripts[code])
        return result


def fasta_lines(source, path):
    for n, line in enumerate(source.file, 1):
        if line.translate(None, _FASTA_BYTES):
            raise Unsupported('%s line %d: a carriage return, control or non-ASCII byte' % (path, n))
        yield line.decode('ascii')


def fasta_records(lines):
    """(tag, accession, description, sequence) of each record, exactly as modules.common.read_fasta reads them."""
    metadata = next(lines, '')
    while metadata != '':
        parts = []
        line = next(lines, '')
        while not line.startswith('>') and not line == '':
            parts.append(line[:-1])
            line = next(lines, '')
        tag = ''
        accession = ''
        description = ''
        if '|' in metadata:
            metadata_parsed = metadata[1:].split('|')
            if 'generic' in metadata_parsed[0]:
                tag = metadata_parsed[0]
            else:
                tag = 'generic_' + metadata_parsed[0]
            if len(metadata_parsed) == 2:
                if ' ' in metadata_parsed[1]:
                    accession = metadata_parsed[1].split(' ')[0]
                    description = metadata_parsed[1].split(' ', 1)[1]
                else:
                    accession = metadata_parsed[1]
            elif len(metadata_parsed) == 3:
                accession = metadata_parsed[1]
                description = metadata_parsed[2]
        else:
            accession = metadata[1:].split(' ')[0]
            if ' ' in metadata:
                description = metadata.split(' ', 1)[1]
        yield tag, accession, description.replace('\n', ''), ''.join(parts)
        metadata = line


class GeneNames:
    """The gene names (or gene ID) the reference collects for a transcript, cached per transcript."""

    def __init__(self, annot):
        self.annot = annot
        self.cache = {}

    def __call__(self, transcript):
        names = self.cache.get(transcript)
        if names is None:
            try:
                feature_attr = self.annot[transcript].attributes
            except Exception as e:      # FeatureNotFoundError, or a key sqlite cannot use: the reference stops
                raise Unsupported('transcript %r: %s: %s' % (transcript, type(e).__name__, e))
            names = ()
            if ('gene_id' in feature_attr):
                names = tuple(feature_attr['gene_name'] if ('gene_name' in feature_attr) else feature_attr['gene_id'])
            self.cache[transcript] = names
        return names


def simplify_block(block, table, genes, fasta_out, header_out, first):
    """Writes one block of records to both outputs."""
    parsed = []
    haplotypes = []
    for tag, accession, description, sequence in block:
        try:
            proteinIDs = description.split('matching_proteins:', 1)[1].split(maxsplit=1)[0]
            proteinPos = description.split('position_within_protein:', 1)[1].split(maxsplit=1)[0]
            proteinStart = description.split('start:', 1)[1].split(maxsplit=1)[0]
            proteinRF = description.split('reading_frame:', 1)[1].split(maxsplit=1)[0]
        except IndexError:
            raise Unsupported('record %s: a header without the fields the reference needs' % accession)
        ids = re.split(r'[;,]', proteinIDs)
        if tag == 'generic_enshap':
            n_before = len(haplotypes)
            haplotypes.extend(protID for protID in ids if protID.startswith(args.haplo_prefix))
            parsed.append((tag, accession, sequence, proteinIDs, proteinPos, proteinStart, proteinRF, ids,
                           len(haplotypes) - n_before))
        else:
            parsed.append((tag, accession, sequence, proteinIDs, proteinPos, proteinStart, proteinRF, ids, 0))
    if haplotypes:
        if table is None:
            raise Unsupported('haplotype records but no haplotype table')
        transcripts = iter(table.transcripts_of(haplotypes))

    fasta = []
    rows = []
    for tag, accession, sequence, proteinIDs, proteinPos, proteinStart, proteinRF, ids, n_haplotypes in parsed:
        if (tag == 'generic_cont'):
            contIDs = [protID for protID in ids if not (protID.startswith('ENST') or protID.startswith(args.haplo_prefix) or protID.startswith(args.var_prefix))]
            fasta.append('>' + accession + ' CONTAMINANT GN=' + ';'.join(contIDs) + '\n' + sequence + '\n')
            rows.append([accession, tag, proteinIDs, proteinPos, proteinStart, proteinRF])
            continue

        if (tag == 'generic_ensref'):
            trIDs = [protID for protID in ids if (protID.startswith('ENST'))]
        elif (tag == 'generic_enshap'):
            trIDs = [next(transcripts) for _ in range(n_haplotypes)]
        elif (tag == 'generic_var'):
            try:
                trIDs = ['ENST' + protID.split('ENST', 1)[1].split('_', 1)[0] for protID in ids if (protID.startswith(args.var_prefix))]
            except IndexError:
                raise Unsupported('record %s: a variant ID without ENST' % accession)
        else:
            trIDs = []

        gIDs = []
        for trID in trIDs:
            gIDs.extend(genes(trID))
        gIDs = list(dict.fromkeys(gIDs))

        if (len(gIDs) > 1):
            for i, gID in enumerate(gIDs):
                fasta.append('>' + accession + '.' + str(i) + ' GN=' + gID + '\n' + sequence + '\n')
                rows.append([accession + '.' + str(i), tag, proteinIDs, proteinPos, proteinStart, proteinRF])
        elif (len(gIDs) == 1):
            fasta.append('>' + accession + ' GN=' + gIDs[0] + '\n' + sequence + '\n')
            rows.append([accession, tag, proteinIDs, proteinPos, proteinStart, proteinRF])
        else:
            fasta.append('>' + accession + ' GN=missing\n' + sequence + '\n')
            rows.append([accession, tag, proteinIDs, proteinPos, proteinStart, proteinRF])

    try:
        # record by record: joining a block into one string makes macOS's allocator keep every such string's pages
        for text in fasta:
            fasta_out.write(text.encode('ascii'))
        # Every value is a str, so to_csv writes (and quotes) each row on its own: blocks join into the whole table.
        header_out.write(pd.DataFrame(data=rows, columns=HEADER_COLUMNS).to_csv(sep='\t', index=False, header=first)
                         .encode('ascii'))
    except UnicodeEncodeError:
        raise Unsupported('a gene name with non-ASCII characters')


def simplify():
    for out in (args.output_file, args.output_file_header):
        for inp in [args.input_file] + ([args.haplo_table] if args.haplo_table else []):
            if os.path.exists(out) and os.path.exists(inp) and os.path.samefile(out, inp):
                raise Unsupported('%s is both an input and an output' % out)
    check_compression(args.output_file_header)
    table = HaplotypeTable(args.haplo_table) if (args.haplo_table) else None
    try:
        annot = gffutils.FeatureDB(args.annotation_db)
    except Exception as e:
        raise Unsupported('annotation DB: %s: %s' % (type(e).__name__, e))
    genes = GeneNames(annot)
    source = Input(args.input_file)
    fasta_out = GzipWriter(args.output_file) if args.output_file.endswith('.gz') else open(args.output_file, 'wb', buffering=1 << 20)
    header_out = Output(args.output_file_header)
    accessions = array('q')
    try:
        records = fasta_records(fasta_lines(source, args.input_file))
        first = True
        while True:
            block = list(islice(records, BLOCK))
            if not block and not first:
                break
            accessions.extend(hash(record[1]) for record in block)
            simplify_block(block, table, genes, fasta_out, header_out, first)
            first = False
        source.close()
    except Unsupported:
        source.kill()
        raise
    finally:
        fasta_out.close()
        header_out.close()
    # The reference keeps one record per accession (the first one's place, the last one's content).
    hashes = np.sort(np.frombuffer(accessions, dtype=np.int64))
    if len(hashes) > 1 and bool((hashes[1:] == hashes[:-1]).any()):
        raise Unsupported('two records with the same accession (or the same hash of it)')


try:
    simplify()
except Unsupported as e:
    to_legacy(__file__, str(e))
