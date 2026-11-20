import argparse
import glob
import os
import subprocess
import sys

import numpy as np

from modules.fallback import ERROR_STATUS, to_legacy
from modules.gzip_writer import GzipWriter

# ProHap-fast (IMPLEMENTATION_NOTES section 2): the same arguments, output and printed lines as the reference
# (legacy/src/vcf_filter_fix.py), with the genotypes handled as bytes by numpy instead of one Python string per
# sample. A line the fast path does not cover (FORMAT other than GT, genotypes other than one character per allele,
# whitespace other than single tabs, ten or more ALT alleles) goes through the reference's own code for one line,
# copied below. Any carriage return or non-ASCII byte hands the whole run to the reference, whose text-mode reading
# treats those differently from bytes. The output is gzip level 1, compressed by a thread (the reference uses level
# 9; the content is the same).
#
# `-o -` (fast only, used by prohap.py's -filter_glob, F3): the filtered VCF goes to stdout uncompressed and the
# messages to stderr; input the fast path does not handle stops the script with exit status 3 instead of handing it to
# the reference, and prohap.py then runs the steps separately.

parser = argparse.ArgumentParser(
        description='Reads the VCF file, parses multi-allelic variants into multiple lines, and filters out variants under the MAF threshold.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input VCF")

parser.add_argument("-af_field", dest="af_field", required=False, type=str,
                    help="Allele Frequency (AF) field name - default \"AF\"", default="AF")

parser.add_argument("-af", dest="min_af", required=False, type=float,
                    help="Allele Frequency (AF) lower threshold - default 0", default=0)

parser.add_argument("-chr", dest="chromosome", required=True,
                    help="chromosome being processed (e.g., 1, 12 or X)")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output VCF; - for uncompressed text on stdout (ProHap-fast only)")

args = parser.parse_args()


def get_MAF(info):
    if ';' + args.af_field + '=' in info:
        return info.split(';' + args.af_field + '=')[1].split(';')[0]
    elif '\t' + args.af_field + '=' in info:
        return info.split('\t' + args.af_field + '=')[1].split(';')[0]
    elif info.startswith(args.af_field + '='):
        return info.split(args.af_field + '=')[1].split(';')[0]

    return "-1"


def reference_line(line, outfile):
    """The reference's loop body for one data line (str), verbatim but for the two counters it returns."""
    total_VCF_entries = 0
    valid_VCF_entries = 0

    ALT = line.split(maxsplit=5)[4]
    INFO = line.split(maxsplit=8)[7]
    GT_IDX = line.split(maxsplit=9)[8].split(':').index('GT')

    # If not thresholding by MAF, ignore any possible formatting errors
    MAF = "0"
    if (args.min_af > 0):
        MAF = get_MAF(INFO)

    total_VCF_entries += 1

    # ProHap does not work with multi-allelic variants -> make a separate line in the output VCF for each allele
    if (',' in ALT):
        for i,allele in enumerate(ALT.split(',')):
            if (args.min_af > 0):
                allele_maf = float(MAF.split(',')[i])
            else:
                allele_maf = 0

            # keep only the i-th allele as the alternatvie (1), mark all the other alleles in the genotypes as 0
            if (allele_maf >= args.min_af):
                invalid_gts = list(range(1,100))    # assuming that a variant will not have more than 100 possible alleles
                invalid_gts.remove(i+1)

                # format the string containing all the genotypes
                GTs = [ GT.split(':')[GT_IDX] for GT in line.split()[9:] ]

                # on sex chromosomes for male individuals, there might be just a single allele and the genotype won't be formated as "x|x" -> keep the valid allele as the first, and add a 0 to the genotype for consistency
                GTs_str = '\t'.join([ ((GT + '|0') if (('|' not in GT) and ('/' not in GT)) else GT) for GT in GTs ])

                # mark all the other alleles as 0
                for gt_id in invalid_gts:
                    GTs_str = GTs_str.replace(str(gt_id) + '|', '0|')
                    GTs_str = GTs_str.replace('|' + str(gt_id), '|0')

                # mark this allele as 1
                GTs_str = GTs_str.replace(str(i+1) + '|', '1|')
                GTs_str = GTs_str.replace('|' + str(i+1), '|1')

                outfile.write(('\t'.join(line.split(maxsplit=4)[:-1] + [allele] + line.split(maxsplit=7)[5:-1])).encode('ascii'))
                outfile.write(('\tMAF=' + str(allele_maf) + '\tGT\t' + GTs_str + '\n').encode('ascii'))
                valid_VCF_entries += 1

    else:
        allele_maf = float(MAF)
        if (allele_maf >= args.min_af):
            # format the string containing all the genotypes
            GTs = [ GT.split(':')[GT_IDX] for GT in line.split()[9:] ]

            # on sex chromosomes for male individuals, there might be just a single allele and the genotype won't be formated as "x|x"
            # -> keep the valid allele as the first, and add a 0 to the genotype for consistency
            GTs_str = '\t'.join([ ((GT + '|0') if (('|' not in GT) and ('/' not in GT)) else GT) for GT in GTs ])

            outfile.write(('\t'.join(line.split(maxsplit=8)[:-1]) + '\tGT\t' + GTs_str + '\n').encode('ascii'))
            valid_VCF_entries +=1

    return total_VCF_entries, valid_VCF_entries


# ---------------------------------------------------------------------------------------------------------------
# genotypes as bytes

TAB, NL, PIPE, SLASH, ZERO = (ord(c) for c in '\t\n|/0')
REF_WORD = np.frombuffer(b'0|0\t', dtype='<u4')[0]
IS_ALLELE = np.zeros(256, dtype=bool)               # one character per allele: a digit or '.'
IS_ALLELE[[ord(c) for c in '0123456789.']] = True
IS_SEP = np.zeros(256, dtype=bool)
IS_SEP[[PIPE, SLASH]] = True
# REMAP[j]: the reference's str.replace chain for ALT allele j (1-9) on one allele of a '|' genotype: j becomes 1, every
# other non-zero digit 0; '0' and '.' stay
REMAP = np.tile(np.arange(256, dtype=np.uint8), (10, 1))
for _j in range(1, 10):
    for _d in range(1, 10):
        REMAP[_j, ord(str(_d))] = ord('1') if _d == _j else ZERO


def normalised_genotypes(line, start):
    """The genotypes of a line (from byte `start` to the end, which is a newline) as the reference writes them: tokens
    of one character per allele, `a` (haploid) written `a|0`. As bytes, ending in a newline; None when a token is
    anything else, or a separator anything but a single tab."""
    size = len(line) - start
    if size % 4 == 0:
        # all diploid: 'a|b\t' words, the last 'a|b\n'
        words = np.frombuffer(line, dtype='<u4', offset=start)
        other = np.flatnonzero(words != REF_WORD)       # the last word always (newline)
        w = words[other]
        end = np.where(other == len(words) - 1, NL, TAB)
        if (IS_ALLELE[w & 0xff] & IS_SEP[(w >> 8) & 0xff] & IS_ALLELE[(w >> 16) & 0xff] & ((w >> 24) == end)).all():
            return line[start:]
    # haploid and diploid tokens mixed (chromosome X)
    a = np.frombuffer(line, dtype=np.uint8, offset=start)
    ends = np.flatnonzero((a == TAB) | (a == NL))
    if len(ends) == 0 or ends[-1] != len(a) - 1:
        return None
    starts = np.empty_like(ends)
    starts[0] = 0
    starts[1:] = ends[:-1] + 1
    lengths = ends - starts
    diploid = lengths == 3
    if not (diploid | (lengths == 1)).all():
        return None
    out = np.empty((len(ends), 4), dtype=np.uint8)
    out[:, 0] = a[starts]
    out[:, 1] = np.where(diploid, a[np.minimum(starts + 1, len(a) - 1)], PIPE)
    out[:, 2] = np.where(diploid, a[np.minimum(starts + 2, len(a) - 1)], ZERO)
    out[:, 3] = TAB
    out[-1, 3] = NL
    if not (IS_ALLELE[out[:, 0]] & IS_SEP[out[:, 1]] & IS_ALLELE[out[:, 2]]).all():
        return None
    return out.tobytes()


def remapped_genotypes(genotypes, allele):
    """Genotypes for ALT allele `allele` (1-9) of a multi-allelic line: only tokens with a '|' are changed."""
    out = np.frombuffer(genotypes, dtype=np.uint8).reshape(-1, 4).copy()
    phased = out[:, 1] == PIPE
    remap = REMAP[allele]
    out[:, 0] = np.where(phased, remap[out[:, 0]], out[:, 0])
    out[:, 2] = np.where(phased, remap[out[:, 2]], out[:, 2])
    return out.tobytes()


def fast_line(line, outfile):
    """One data line (bytes, ending in a newline) as the reference handles it; None (nothing written) when the line
    needs the reference's own code."""
    tabs = []
    t = -1
    for _ in range(9):
        t = line.find(b'\t', t + 1)
        if t < 0:
            return None
        tabs.append(t)
    head = line[:tabs[8]].decode('ascii')
    # the reference splits on any whitespace: the fast path needs nine non-empty fields between single tabs
    if len(head.split()) != 9 or head.count('\t') != 8:
        return None
    fields = head.split('\t')
    if fields[8] != 'GT':
        return None
    ALT = fields[4]
    INFO = fields[7]
    MAF = "0"
    if (args.min_af > 0):
        MAF = get_MAF(INFO)

    if ',' not in ALT:
        allele_maf = float(MAF)
        if allele_maf >= args.min_af:
            genotypes = normalised_genotypes(line, tabs[8] + 1)
            if genotypes is None:
                return None
            if len(genotypes) == len(line) - tabs[8] - 1:
                outfile.write(line)                    # FORMAT is GT and nothing to normalise: the line as it is
            else:
                outfile.write(line[:tabs[7]] + b'\tGT\t' + genotypes)
            return 1, 1
        return 1, 0

    alleles = ALT.split(',')
    if len(alleles) > 9:
        return None
    lines = []
    genotypes = None
    for i, allele in enumerate(alleles):
        if (args.min_af > 0):
            allele_maf = float(MAF.split(',')[i])
        else:
            allele_maf = 0
        if allele_maf >= args.min_af:
            if genotypes is None:
                genotypes = normalised_genotypes(line, tabs[8] + 1)
                if genotypes is None:
                    return None
            lines.append(('\t'.join(fields[:4] + [allele] + fields[5:7]) + '\tMAF=' + str(allele_maf)
                          + '\tGT\t').encode('ascii') + remapped_genotypes(genotypes, i + 1))
    for out in lines:
        outfile.write(out)
    return 1, len(lines)


class StdoutWriter:
    """The filtered VCF as text on stdout, for prohap.py to read through a pipe."""

    def __init__(self):
        self.out = sys.stdout.buffer
        self.pending = []
        self.pending_size = 0

    def write(self, data):
        self.pending.append(data)
        self.pending_size += len(data)
        if self.pending_size >= 1 << 20:
            self.out.write(b''.join(self.pending))
            self.pending = []
            self.pending_size = 0

    def close(self):
        self.out.write(b''.join(self.pending))
        self.pending = []
        self.out.flush()


def unsupported(reason):
    for proc in processes:
        if proc.poll() is None:
            proc.kill()
    if to_stdout:
        sys.stdout.flush()
        print('%s (%s): %s' % (STOP_MESSAGE, 'vcf_filter_fix.py', reason), file=sys.stderr, flush=True)
        os._exit(ERROR_STATUS)
    to_legacy(__file__, reason)


def open_input(path):
    if path.endswith('.gz'):
        proc = subprocess.Popen(['gzip', '-dc', path], stdout=subprocess.PIPE, bufsize=1 << 20)
        processes.append(proc)
        return proc.stdout, proc
    return open(path, 'rb', buffering=1 << 20), None


processes = []

# if there are multiple files for the same chromosome (specified, e.g., as "chr22_*.vcf"), process all of them in the alphabetical order
input_files = sorted(glob.glob(args.input_file))

# the output, gzipped as it is written (or text on stdout)
to_stdout = args.output_file == '-'
STOP_MESSAGE = 'ProHap-fast: the filter stops on input it does not handle'
outfile = StdoutWriter() if to_stdout else GzipWriter(args.output_file)

# Remember the number of original entries and the number of filtered ones
total_VCF_entries = 0
valid_VCF_entries = 0
messages = []

for i in range(len(input_files)):
    vcf_file, reader = open_input(input_files[i])
    line_number = 0

    # the header: every line before the first data line that starts with '#'; written for the first file only
    line = vcf_file.readline()
    while line != b'' and line.startswith(b'#'):
        line_number += 1
        if b'\r' in line or not line.isascii():
            unsupported('%s line %d: carriage return or non-ASCII byte' % (input_files[i], line_number))
        if i == 0:
            outfile.write(line)
        line = vcf_file.readline()

    if line == b'':
        messages.append(input_files[i] + ' is empty!')

    while line != b'':
        line_number += 1
        if b'\r' in line or not line.isascii():
            unsupported('%s line %d: carriage return or non-ASCII byte' % (input_files[i], line_number))
        if not line.endswith(b'\n'):
            line += b'\n'
        counts = fast_line(line, outfile)
        if counts is None:
            counts = reference_line(line.decode('ascii'), outfile)
        total_VCF_entries += counts[0]
        valid_VCF_entries += counts[1]
        line = vcf_file.readline()

    vcf_file.close()
    if reader is not None and reader.wait() != 0:
        unsupported('gzip could not decompress ' + input_files[i])

outfile.close()

log = sys.stderr if to_stdout else sys.stdout
for message in messages:
    print(message, file=log)
print(('Chr ' + args.chromosome + ':'), 'Original VCF(s):', total_VCF_entries, 'lines', file=log)
print(('Chr ' + args.chromosome + ':'), 'Filtered VCF:', valid_VCF_entries, 'lines', file=log)
