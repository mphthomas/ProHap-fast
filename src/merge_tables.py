import argparse

from modules.fallback import Unsupported, to_legacy
from modules.tsv_stream import Output, Table, check_compression, combined_dtype

# ProHap-fast (IMPLEMENTATION_NOTES section 4): the same arguments and output as the reference
# (legacy/src/merge_tables.py), which reads every table into pandas, joins them with pd.concat and writes the result
# with to_csv. This script streams the tables chunk by chunk instead (modules/tsv_stream.py), so it holds one chunk
# of rows rather than all of them. Tables it does not handle (see tsv_stream), different headers, or columns whose
# chunks pandas would join by rewriting values go to the reference. A gzipped output is written at level 1 (the
# reference uses pandas' default, 9; the content is the same).

parser = argparse.ArgumentParser(
    description='Merge all the haplotype tables into one TSV file.')

parser.add_argument("-i", dest="input_filenames", required=True,
                    help="input files, comma-separated list")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output TSV file")

args = parser.parse_args()


def merge(inputs, output_file):
    check_compression(output_file)
    names = None
    output = None
    merged = None       # per column: the dtypes the tables give it
    table = None
    try:
        for path in inputs:
            table = Table(path)
            if names is None:
                names = table.names
                merged = [set() for _ in names]
                output = Output(output_file)
                output.write(('\t'.join(names) + '\n').encode('ascii'))
            elif table.names != names:
                raise Unsupported(path + ': the columns differ from those of ' + inputs[0])
            table.read(output.write_rows)
            if table.rows > 0:
                for i, name in enumerate(names):
                    merged[i].add(combined_dtype(table.dtypes[i], name))
            table = None
        for i, name in enumerate(names):
            combined_dtype(merged[i], name)
    except Unsupported:
        if table is not None:
            table.kill()
        raise
    finally:
        if output is not None:
            output.close()


try:
    merge(args.input_filenames.split(","), args.output_file)
except Unsupported as e:
    to_legacy(__file__, str(e))
