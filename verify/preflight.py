"""Preflight checks of the phased VCF chunks that verify is about to use.

Everything here only warns: verify still runs both arms. Each chunk's header is read in full, and its records only
up to MAX_RECORDS or MAX_BYTES (a record of a 500k-sample VCF is about 2 MB), so the checks of records and of chunk
order are a sample of the file, not a proof.
"""
import collections
import gzip
import os
import re
import shutil
import subprocess

import vcfcut

MAX_RECORDS = 1000
MAX_BYTES = 256 << 20

# one character per allele: 0|1, .|0, haploid 1
_MULTI_CHAR = re.compile(rb'[0-9.][0-9.]')
_OTHER_CHAR = re.compile(rb'[^0-9.|/\t\n]')
_NUMERIC = re.compile(r'[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?')


class Chunk:
    """What one scan of a VCF chunk found."""

    def __init__(self, path):
        self.path = path
        self.samples = None         # sample names from the #CHROM line, None without one
        self.records = 0            # records scanned
        self.eof = False            # every record was scanned
        self.chroms = []            # CHROM values in order of first appearance
        self.first_pos = None
        self.last_pos = None        # last scanned position (the file's last one only when eof)
        self.unsorted = 0
        self.gt_not_first = 0       # FORMAT with GT, but not as the first field
        self.no_gt = 0
        self.unphased = 0           # a genotype with '/'
        self.multi_char = 0         # an allele written with more than one character
        self.other_char = 0         # any other byte in the genotype columns (CR, letters, ...)
        self.dot_ids = 0


def open_vcf(path):
    """Open a chunk as the legacy filter does: gzip when the name ends in .gz, else plain text."""
    return gzip.open(path, 'rb') if path.endswith('.gz') else open(path, 'rb')


def names(line):
    """The sample names of a #CHROM line."""
    return [n.decode('utf-8', 'replace') for n in vcfcut.sample_names(line)]


def read_samples(path):
    """Sample names of a chunk's #CHROM line; None without a #CHROM line."""
    with open_vcf(path) as fh:
        for line in fh:
            if not line.startswith(b'#'):
                return None
            if line.startswith(b'#CHROM'):
                return names(line)
    return None


def scan(path, max_records=None, max_bytes=None):
    max_records = MAX_RECORDS if max_records is None else max_records
    max_bytes = MAX_BYTES if max_bytes is None else max_bytes
    chunk = Chunk(path)
    read = 0
    with open_vcf(path) as fh:
        for line in fh:
            if line.startswith(b'#'):
                if line.startswith(b'#CHROM'):
                    chunk.samples = names(line)
                continue
            if chunk.records >= max_records or read >= max_bytes:
                return chunk
            read += len(line)
            chunk.records += 1
            fields = line.split(b'\t', 9)
            if len(fields) < 2:
                continue
            chrom = fields[0].decode('utf-8', 'replace')
            if chrom not in chunk.chroms:
                chunk.chroms.append(chrom)
            try:
                pos = int(fields[1])
            except ValueError:
                pos = None
            if pos is not None:
                if chunk.first_pos is None:
                    chunk.first_pos = pos
                if chunk.last_pos is not None and pos < chunk.last_pos:
                    chunk.unsorted += 1
                chunk.last_pos = pos
            if len(fields) > 2 and fields[2] == b'.':
                chunk.dot_ids += 1
            if len(fields) < 10:
                continue
            fmt = fields[8].split(b':')
            if b'GT' not in fmt:
                chunk.no_gt += 1
                continue
            if fmt[0] != b'GT':
                chunk.gt_not_first += 1
                continue
            gts = fields[9]
            if len(fmt) > 1:            # GT:DS, ...: the filter keeps the GT subfield of each token
                gts = b'\t'.join(token.split(b':', 1)[0] for token in gts.split(b'\t'))
            if b'/' in gts:
                chunk.unphased += 1
            if _MULTI_CHAR.search(gts):
                chunk.multi_char += 1
            if _OTHER_CHAR.search(gts):
                chunk.other_char += 1
    chunk.eof = True
    return chunk


def records_after(path, chrom, pos, tabix):
    """Whether an indexed chunk has a record at or after chrom:pos; None when it has no index."""
    if not tabix or not any(os.path.exists(path + ext) for ext in ('.tbi', '.csi')):
        return None
    proc = subprocess.run([tabix, path, chrom + ':' + str(pos)], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    if proc.returncode != 0:
        return None
    return len(proc.stdout) > 0


def _examples(items, n=5):
    items = list(items)
    return ', '.join(items[:n]) + (', ...' if len(items) > n else '')


def check(chunks_by_chrom, metadata_names, tabix=None):
    """Warnings and notes for the chunks of each chromosome. metadata_names: the set of sample names the legacy code
    reads from sample_metadata_file. Returns (warnings, notes), lists of strings."""
    tabix = tabix or shutil.which('tabix')
    warnings, notes = [], []
    seen_numeric, seen_missing = set(), set()
    dot_ids = []
    for chrom, paths in chunks_by_chrom.items():
        if not paths:
            warnings.append('chromosome %s: no VCF matches the input pattern, so it gets no haplotypes' % chrom)
            continue
        bcf = [p for p in paths if p.endswith('.bcf')]
        if bcf:
            notes.append('%s: BCF chunks are not scanned (restrict_to_exons reads them with bcftools)'
                         % ', '.join(os.path.relpath(p) for p in bcf))
            paths = [p for p in paths if not p.endswith('.bcf')]
            if not paths:
                continue
        scanned = [scan(p) for p in paths]
        first = scanned[0]
        for c in scanned:
            name = os.path.relpath(c.path)
            sampled = '' if c.eof else ' (first %d records)' % c.records
            if c.samples is None:
                warnings.append('%s has no #CHROM header line' % name)
                continue
            if c is not first and c.samples != first.samples:
                warnings.append('%s lists other samples than %s: the legacy filter keeps only the first chunk\'s '
                                'header, so the genotype columns no longer match the sample names'
                                % (name, os.path.relpath(first.path)))
            dupes = sorted(s for s, n in collections.Counter(c.samples).items() if n > 1)
            if dupes:
                warnings.append('%s lists samples more than once: %s' % (name, _examples(dupes)))
            if c.records == 0:
                notes.append('%s has no records' % name)
                continue
            if c.no_gt:
                warnings.append('%s: %d records%s have no GT in FORMAT' % (name, c.no_gt, sampled))
            if c.gt_not_first:
                warnings.append('%s: %d records%s have FORMAT fields before GT. ProHap reads only GT; the fast filter '
                                'runs the legacy code on each such record, at the legacy speed. Put GT first, or keep '
                                'only GT (bcftools annotate --no-version -x ^FORMAT/GT; this removes header lines, '
                                'which renumbers records with ID ".")' % (name, c.gt_not_first, sampled))
            if c.unphased:
                warnings.append('%s: %d records%s have unphased genotypes (/). ProHap needs phased genotypes (|)'
                                % (name, c.unphased, sampled))
            if c.multi_char:
                warnings.append('%s: %d records%s have alleles written with more than one character (e.g. 10|0); '
                                'the fast implementation runs the legacy code on such files'
                                % (name, c.multi_char, sampled))
            if c.other_char:
                warnings.append('%s: %d records%s have genotype columns with characters other than 0-9 . | / '
                                '(e.g. CR line ends)' % (name, c.other_char, sampled))
            if c.unsorted:
                warnings.append('%s is not sorted by position: %d records%s come before the one above them'
                                % (name, c.unsorted, sampled))
            if len(c.chroms) > 1:
                warnings.append('%s holds more than one chromosome: %s' % (name, ', '.join(c.chroms)))
            if c.dot_ids:
                dot_ids.append('%s (%d of %d%s)' % (name, c.dot_ids, c.records, '' if c.eof else ' scanned'))
            numeric = [s for s in c.samples if _NUMERIC.fullmatch(s)]
            if numeric and not seen_numeric.issuperset(numeric):
                seen_numeric.update(numeric)
                warnings.append('%s: %d sample names are numbers (%s). pandas reads such names in the metadata as '
                                'numbers, and the legacy code then drops those samples'
                                % (name, len(numeric), _examples(numeric)))
            if metadata_names is not None:
                missing = [s for s in c.samples if s not in metadata_names]
                if missing and not seen_missing.issuperset(missing):
                    seen_missing.update(missing)
                    warnings.append('%s: %d of %d samples are not in sample_metadata_file (%s); ProHap leaves them '
                                    'out and lists them in log/prohap_chr*.log'
                                    % (name, len(missing), len(c.samples), _examples(missing)))
        # chunk order: the legacy filter reads the chunks in the alphabetical order of their names
        for a, b in zip(scanned, scanned[1:]):
            if a.first_pos is None or b.first_pos is None:
                continue
            na, nb = os.path.relpath(a.path), os.path.relpath(b.path)
            if b.first_pos < a.first_pos or (a.chroms and b.chroms and a.chroms[0] != b.chroms[0]):
                warnings.append('%s sorts before %s but starts later (position %d vs %d); chunk names must sort in '
                                'genomic order' % (na, nb, a.first_pos, b.first_pos))
            elif a.eof:
                if a.last_pos >= b.first_pos:
                    warnings.append('%s overlaps %s: it ends at %d, after %s starts at %d'
                                    % (na, nb, a.last_pos, nb, b.first_pos))
            else:
                overlap = records_after(a.path, a.chroms[-1], b.first_pos, tabix)
                if overlap:
                    warnings.append('%s overlaps %s: it has records at or after position %d, where %s starts'
                                    % (na, nb, b.first_pos, nb))
                elif overlap is None:
                    notes.append('could not check whether %s overlaps %s: index %s with tabix to check'
                                 % (na, nb, na))
    if dot_ids:
        notes.append('records with ID ".": %s. ProHap names them after their line number in the filtered VCF, so '
                     'their names change when records before them are added or removed. verify keeps every record; '
                     'bcftools annotate --set-id makes the names stable' % _examples(dot_ids, 3))
    return warnings, notes
