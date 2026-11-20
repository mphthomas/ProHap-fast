import argparse

from modules.fallback import Unsupported, to_legacy
from modules.tsv_stream import Output, Table, check_compression, combined_dtype

# ProHap-fast (IMPLEMENTATION_NOTES section 4): the same arguments and outputs as the reference
# (legacy/src/haplo_extract_sample_names.py), which reads the whole table into pandas and writes it twice with
# to_csv. This script streams the table chunk by chunk instead (modules/tsv_stream.py) and writes both outputs from
# the same re-formatted fields, so it holds one chunk of rows rather than all of them. Tables it does not handle (see
# tsv_stream), or columns whose chunks pandas would join by rewriting values, go to the reference. Gzipped outputs are
# written at level 1 (the reference uses pandas' default, 9; the content is the same).

parser = argparse.ArgumentParser(description='Reads a TSV file produced by ProHap, extracts the "samples" column to a separate file and removes it from the original')

parser.add_argument("-hap_tsv", dest="haplo_table", required=True,
                    help="input haplotype table", metavar="FILE")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output TSV file", metavar="FILE")

parser.add_argument("-samples", dest="output_file_samples", required=True,
                    help="output file (tab-separated) for the list of samples per haplotype", metavar="FILE")

args = parser.parse_args()


def extract(input_file, output_file, output_file_samples):
    check_compression(output_file)
    check_compression(output_file_samples)
    table = Table(input_file)
    if 'HaplotypeID' not in table.names:
        table.kill()
        raise Unsupported(input_file + ': no HaplotypeID column')
    column = table.names.index('HaplotypeID')
    outputs = []
    try:
        table_out = Output(output_file)
        outputs.append(table_out)
        samples_out = Output(output_file_samples)
        outputs.append(samples_out)
        table_out.write(('\t'.join(table.names[:-1]) + '\n').encode('ascii'))
        samples_out.write(b'HaplotypeID\tsamples\n')

        def write(prefixes, samples):
            table_out.write(b'\n'.join(prefixes) + b'\n')
            samples_out.write_rows([p.split(b'\t', column + 1)[column] for p in prefixes], samples)
        table.read(write)
        for dtypes, name in zip(table.dtypes, table.names):
            combined_dtype(dtypes, name)
    except Unsupported:
        table.kill()
        raise
    finally:
        for output in outputs:
            output.close()


try:
    extract(args.haplo_table, args.output_file, args.output_file_samples)
except Unsupported as e:
    to_legacy(__file__, str(e))
