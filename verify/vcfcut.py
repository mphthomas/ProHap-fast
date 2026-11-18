"""Cut a VCF to some of its samples, changing nothing else.

    python3 vcfcut.py SRC DST SAMPLES_FILE

Writes SRC to DST with its first nine columns and the columns of the samples named in SAMPLES_FILE (one per line),
in SRC's order. Every other byte stays as it is: header lines, records, INFO text, line ends. So records keep their
line numbers, which name the records whose ID is ".". Files are told apart by name, as the legacy filter does: SRC
is read as gzip when its name ends in .gz, and DST is written as gzip when its name does.

bcftools view -s would re-write floating-point INFO values, add header lines the header lacks, and refuse records
with INFO tags the header does not define.
"""
import gzip
import operator
import sys

FIXED = 9   # CHROM POS ID REF ALT QUAL FILTER INFO FORMAT


def split_end(line):
    """(line without its line end, the line end)."""
    if line.endswith(b'\r\n'):
        return line[:-2], b'\r\n'
    if line.endswith(b'\n'):
        return line[:-1], b'\n'
    return line, b''


def sample_names(line):
    """The sample names of a #CHROM line."""
    return split_end(line)[0].split(b'\t')[FIXED:]


def cut(src, dst, keep):
    """Cut src to the samples in keep (a set of bytes). Returns the number of samples kept."""
    pick = None
    with (gzip.open(src, 'rb') if src.endswith('.gz') else open(src, 'rb')) as fin, \
            (gzip.open(dst, 'wb', compresslevel=1) if dst.endswith('.gz') else open(dst, 'wb')) as fout:
        for line in fin:
            if pick is None:
                if line.startswith(b'#CHROM'):
                    names = sample_names(line)
                    columns = [FIXED + i for i, name in enumerate(names) if name in keep]
                    if not columns:
                        raise ValueError(src + ' has none of the selected samples')
                    pick = operator.itemgetter(*range(FIXED), *columns)
                    width = FIXED + len(names)
                elif not line.startswith(b'#'):
                    raise ValueError(src + ' has no #CHROM header line')
                else:
                    fout.write(line)
                    continue
            body, end = split_end(line)
            fields = body.split(b'\t')
            if len(fields) != width:
                raise ValueError('%s: a line has %d columns where the header has %d' % (src, len(fields), width))
            fout.write(b'\t'.join(pick(fields)) + end)
    if pick is None:
        raise ValueError(src + ' has no #CHROM header line')
    return len(columns)


def main(argv):
    if len(argv) != 3:
        sys.exit(__doc__)
    src, dst, samples_file = argv
    with open(samples_file, 'rb') as fh:
        keep = {name for name in fh.read().split(b'\n') if name}
    try:
        cut(src, dst, keep)
    except ValueError as e:
        sys.exit('vcfcut: ' + str(e))


if __name__ == '__main__':
    main(sys.argv[1:])
