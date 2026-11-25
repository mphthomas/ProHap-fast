'''
Main file of the ProHap tool

Creates a database of CDS + protein haplotypes, and a fasta file of protein haplotype sequences.

ProHap-fast (IMPLEMENTATION_NOTES section 3): the same arguments and output as the reference (legacy/src/prohap.py),
without its per-transcript tables of one column per sample. The VCF is read once (modules/sweep.py); each transcript
is handed to a worker as soon as the reference's queue would release it, and the worker computes its haplotypes
(modules/haplotypes.py) and protein haplotypes (modules/proteins.py). This process writes the results in the order of
the transcripts (all_transcripts), numbering haplotypes and FASTA accessions over the whole chromosome as the
reference does. Input the fast path does not reproduce exactly hands the run to the reference (modules/fallback.py).

With -filter_glob (and -filter_af, -filter_af_field) instead of -i, the filter runs fused (F3, modules/fused.py): the
chunks are filtered by vcf_filter_fix.py in a process of its own and read from it through a pipe, so the filtered VCF
is neither written nor read back.
'''

import argparse
import gc
import gzip
import os
import queue
import sys
from multiprocessing import Pool

import gffutils
import pandas as pd
from numpy import int64
from pandas.io.common import infer_compression

from modules.cdna import read_sequences
from modules.fallback import Unsupported, cannot_fall_back, to_legacy
from modules.fused import FilterReader, to_legacy_fused
from modules.legacy_import import legacy
from modules.pandas_emulation import emulate_id_column, pandas_chunk_rows
from modules.sweep import FIXED_COLUMNS, RecordParser, VCFReader, assign_records, exon_index, read_header
from modules.haplotypes import Cohort, transcript_rows
from modules.fixes import selected as selected_fixes
from modules.proteins import Settings, Transcript, process_transcript

check_open_file = legacy.common.check_open_file
empty_output = legacy.process_haplotypes.empty_output
result_columns = legacy.process_haplotypes.result_columns

parser = argparse.ArgumentParser(
        description='Creates a database of CDS + protein haplotypes, and a fasta file of protein haplotype sequences.')

parser.add_argument("-i", dest="input_vcf", required=False,
                    help="input VCF file (the filtered VCF); or -filter_glob", metavar="FILE",
                    type=lambda x: check_open_file(parser, x))

parser.add_argument("-filter_glob", dest="filter_glob", required=False,
                    help="ProHap-fast: the unfiltered VCF chunks (glob), filtered on the fly as vcf_filter_fix.py does; "
                         "replaces -i")

parser.add_argument("-filter_af", dest="filter_af", required=False, type=float, default=0,
                    help="with -filter_glob: vcf_filter_fix.py's -af")

parser.add_argument("-filter_af_field", dest="filter_af_field", required=False, default="AF",
                    help="with -filter_glob: vcf_filter_fix.py's -af_field")

parser.add_argument("-db", dest="annotation_db", required=True,
                    help="DB file created by gffutils from GTF")

parser.add_argument("-s", dest="samples_filename", required=True,
                    help="tab-separated file with sample information (must include 'Sample name' and 'Sex' columns)")

parser.add_argument("-af", dest="min_af", required=False, type=float,
                    help="Allele Frequency (AF) lower threshold - default 0", default=0)

parser.add_argument("-cdna", dest="cdnas_fasta", required=True,
                    help="input cDNA fasta file")

parser.add_argument("-transcripts", dest="transcript_list", required=True,
                    help="list of transcript IDs, provided in a CSV file", metavar="FILE")

parser.add_argument("-require_start", dest="require_start", required=False, type=int,
                    help="flag: require annotation of the start codon, set to 0 to disable; default: 1", default=1)

parser.add_argument("-ignore_UTR", dest="ignore_UTR", required=False, type=int,
                    help="flag (0 or 1): ignore variation in the UTR sequences, do not add UTR translation to proteins; default: 1", default=1)

parser.add_argument("-skip_start_lost", dest="skip_start_lost", required=False, type=int,
                    help="flag (0 or 1): ignore haplotypes where the start codon is lost; default: 1", default=1)

parser.add_argument("-force_rf", dest="force_rf", required=False,
                    help="Force the most likely reading frame when start codon is not annotated or lost due to mutation, set to 0 to disable; default: 1", default=1)

parser.add_argument("-x_par1_to", dest="x_par1_to", required=False, type=int64,
                    help="end location of the 1st pseudoautosomal region on chromosome X; default: 2,781,479", default=2781479)

parser.add_argument("-x_par2_from", dest="x_par2_from", required=False, type=int64,
                    help="start location of the 2nd pseudoautosomal region on chromosome X; default: 155,701,383", default=155701383)

parser.add_argument("-chr", dest="chromosome", required=True,
                    help="chromosome being processed (e.g., 1, 12 or X)")

parser.add_argument("-threads", dest="threads", required=False, type=int,
                    help="number of threads to use; default: 4", default=4)

parser.add_argument("-min_hap_freq", dest="min_freq", required=False, type=float,
                    help="Minimum frequency of a haplotype to be reported in the result (specify -1 to use count threshold instead); default: -1", default=-1)

parser.add_argument("-min_hap_count", dest="min_hap_count", required=False, type=int,
                    help="Minimum count of occurrences of a haplotype to be reported in the result (used only if frequency threshold is not specified); default: 0", default=0)

parser.add_argument("-tag", dest="fasta_tag", required=False,
                    help="tag for FASTa file entries", default='generic_enshap')

parser.add_argument("-id_prefix", dest="haplo_id_prefix", required=False,
                    help="prefix for the haplotype identifier; default: haplo_", default='haplo_')

parser.add_argument("-acc_prefix", dest="accession_prefix", required=False,
                    help="prefix for FASTA file entries accession; default: enshap", default='enshap')

parser.add_argument("-log", dest="log_file", required=False,
                    help="output log file; default: prohap.log", default="prohap.log")

parser.add_argument("-tmp_dir", dest="tmp_dir", required=False,
                    help="directory for temporary files; default: tmp", default="tmp")

parser.add_argument("-output_csv", dest="output_file", required=True,
                    help="output CSV file")

parser.add_argument("-output_fasta", dest="output_fasta", required=True,
                    help="output FASTA file")

parser.add_argument("-output_cdna_fasta", dest="output_cdna_fasta", required=False, default="",
                    help="output cDNA FASTA file (optional; default: none)")

parser.add_argument("-fix_bugs", dest="fix_bugs", required=False, type=int, choices=(0, 1), default=1,
                    help="ProHap-fast: 1 applies the fixes of the reference's haplotype bugs (IMPLEMENTATION_NOTES "
                         "section 9), 0 writes the reference's output byte for byte; default: 1")

TSV_BATCH_ROWS = 20000          # rows per DataFrame.to_csv call; every column has one type, whatever the batch
IN_FLIGHT_PER_WORKER = 4        # transcripts handed to the workers and not yet back, per worker

# the state the workers need, set before the Pool forks them
TRANSCRIPTS = []
COHORT = None
SETTINGS = None
CHUNK_ROWS = None


def without_fix_bugs(argv):
    """argv without -fix_bugs and its value, which the reference does not take."""
    out = []
    skip = False
    for arg in argv:
        if skip:
            skip = False
        elif arg == '-fix_bugs':
            skip = True
        elif not arg.startswith('-fix_bugs='):
            out.append(arg)
    return out


def haplotypes_of_transcript(index, records):
    """Worker: get_haplotypes, then process_transcript_haplotypes, for one transcript with at least one record.
    Returns ('ok', index, number of table rows, number of removed copies, removed samples, result rows, protein
    sequences, cDNA sequences), or ('unsupported', index, reason)."""
    t = TRANSCRIPTS[index]
    try:
        ids = emulate_id_column([r[2] for r in records], CHUNK_ROWS)
        autosomal = COHORT.is_autosomal(t.start, t.end)
        rows, removed = transcript_rows(records, ids, autosomal, COHORT, SETTINGS.fixes)
        table, proteins, cdnas = process_transcript(t, list(enumerate(rows)), SETTINGS, COHORT)
        labels = COHORT.labels
        removed_text = ';'.join([labels[c] for c in removed.tolist()])
        return 'ok', index, len(rows), len(removed), removed_text, table, proteins, cdnas
    except Unsupported as e:
        return 'unsupported', index, str(e)


def transcript_info(annotations_db, feature, cdnas):
    exons = [exon for exon in annotations_db.children(feature, featuretype='exon', order_by='start')]
    start_codons = [sc for sc in annotations_db.children(feature, featuretype='start_codon', order_by='start')]
    stop_codons = [sc for sc in annotations_db.children(feature, featuretype='stop_codon', order_by='start')]
    try:
        biotype = feature['transcript_biotype'][0]
    except:     # noqa: E722  (the reference's bare except)
        biotype = '-'
    exon_spans = [_Span(e.start, e.end) for e in exons]
    t = Transcript(feature.id, feature.start, feature.end, feature.strand, exon_spans,
                   start_codons[0].start if start_codons else None, stop_codons[0].start if stop_codons else None,
                   biotype, cdnas.get(feature.id))
    t.exon_starts, t.exon_max_ends = exon_index([(e.start, e.end) for e in exons])
    return t


class _Span:
    __slots__ = ('start', 'end')

    def __init__(self, start, end):
        self.start = start
        self.end = end


class Outputs:
    """The output files, written as results arrive in transcript order."""

    def __init__(self, args, prefix):
        self.args = args
        self.prefix = prefix
        self.created = []
        compression = infer_compression(args.output_file, 'infer')
        if compression not in (None, 'gzip'):
            raise Unsupported('output table compressed as ' + compression)
        self.tsv = self._open(args.output_file, compression == 'gzip')
        self.tsv.write('\t'.join(result_columns) + '\n')
        self.fasta = self._open(args.output_fasta, args.output_fasta.endswith('.gz'))
        self.cdna = None
        if len(args.output_cdna_fasta) > 0:
            self.cdna = self._open(args.output_cdna_fasta, args.output_cdna_fasta.endswith('.gz'))
        self.log = self._open(args.log_file, False)
        self.log.write('TranscriptID\t#Removed\tRemovedSamples\n')
        self.rows = []
        self.offset = 0         # haplotype rows so far, the reference haplotypes and the rows dropped later included
        self.accession = 0

    def _open(self, path, compressed):
        handle = gzip.open(path, 'wt', compresslevel=1) if compressed else open(path, 'w')
        self.created.append((path, handle))
        return handle

    def write(self, transcript_id, n_rows, n_removed, removed_text, table, proteins, cdnas):
        name = lambda i: self.prefix + '_' + hex(self.offset + i)[2:]       # noqa: E731
        if n_rows > 0:
            self.log.write(transcript_id + '\t' + str(n_removed) + '\t' + removed_text + '\n')
        for fields in table:
            fields[3] = name(fields[3])
            self.rows.append(fields)
        if len(self.rows) >= TSV_BATCH_ROWS:
            self.flush()
        for seq in proteins:
            accession = self.args.accession_prefix + '_' + hex(self.accession)[2:]
            description = 'matching_proteins:' + ';'.join([name(h) for h in seq['haplotypes']]) + ' start:' + \
                str(seq['start']) + ' reading_frame:' + ';'.join(seq['rfs'])
            self.fasta.write('>' + self.args.fasta_tag + '|' + accession + '|' + description + '\n')
            self.fasta.write(seq['sequence'] + '\n')
            self.accession += 1
        if self.cdna is not None:
            for seq in cdnas:
                self.cdna.write('>' + ';'.join([name(h) for h in seq['haplotypes']]) + ' start:' + str(seq['start'])
                                + '\n')
                self.cdna.write(seq['sequence'] + '\n')
        self.offset += n_rows

    def flush(self):
        if self.rows:
            pd.DataFrame(self.rows, columns=result_columns).to_csv(self.tsv, sep='\t', header=False, index=False)
            self.rows = []

    def close(self):
        self.flush()
        for _, handle in self.created:
            handle.close()

    def discard(self):
        for path, handle in self.created:
            try:
                handle.close()
            except Exception:
                pass
            if os.path.exists(path):
                os.remove(path)


def run(args):
    global TRANSCRIPTS, COHORT, SETTINGS, CHUNK_ROWS
    state = {'pool': None, 'reader': None, 'outputs': None}

    fused = args.filter_glob is not None
    try:
        if fused:
            print('[ProHap] Computing protein haplotypes from', args.filter_glob, '(filtered on the fly)')
            input_vcf = None
        else:
            print('[ProHap] Computing protein haplotypes from', args.input_vcf.name)
            input_vcf = args.input_vcf.name
            args.input_vcf.close()
        if args.min_af > 0:
            raise Unsupported('-af above 0 (the workflow does not filter again here)')

        print (('Chr ' + args.chromosome + ':'), 'Reading', args.annotation_db)
        annotations_db = gffutils.FeatureDB(args.annotation_db)

        print (('Chr ' + args.chromosome + ':'), 'Reading', args.transcript_list)
        transcript_df = pd.read_csv(args.transcript_list)
        transcript_df['chromosome'] = transcript_df['chromosome'].apply(lambda x: str(x))
        transcript_list = transcript_df[transcript_df['chromosome'] == args.chromosome]['transcriptID'].tolist()

        print (('Chr ' + args.chromosome + ':'), 'Reading', args.samples_filename)
        samples_df = pd.read_csv(args.samples_filename, sep='\t')

        print (('Chr ' + args.chromosome + ':'), 'Assigning annotations to transcripts.')
        all_transcripts = []
        for transcript_id in transcript_list:
            feature = annotations_db[transcript_id]
            if (args.require_start):
                start_codons = [ sc for sc in annotations_db.children(feature, featuretype='start_codon', order_by='start') ]
                if (len(start_codons) > 0):
                        all_transcripts.append(feature)
            else:
                all_transcripts.append(feature)
        all_transcripts.sort(key=lambda x: x.start)

        # the header of the VCF, and whether it has any record. The filter's pipe can be read only once: it stays
        # open for the sweep.
        reader = FilterReader(args) if fused else VCFReader(input_vcf)
        state['reader'] = reader
        if fused:
            first_line_number, columns, first_line = read_header(reader)
            if first_line == b'':
                reader.finish()
                state['reader'] = None
        else:
            try:
                first_line_number, columns, first_line = read_header(reader)
            finally:
                reader.close()
                state['reader'] = None
        if first_line == b'':
            print (('Chr ' + args.chromosome + ':'), 'Assigning variants to transcripts.')
            print(('Chr ' + args.chromosome + ':'), 'VCF file is empty, creating empty output files.')
            empty_output(args.output_file, args.output_fasta)
            return

        if len(all_transcripts) == 0:
            raise Unsupported('no transcript of chromosome %s to process (the reference stops)' % args.chromosome)
        ids = [f.id for f in all_transcripts]
        if len(set(ids)) != len(ids):
            raise Unsupported('a transcript is listed twice')
        metadata_sample_names = set(samples_df['Sample name'])
        if any(c in metadata_sample_names for c in FIXED_COLUMNS):
            raise Unsupported('the metadata names a sample like a VCF column')
        sample_ids = [sample for sample in columns if (sample in metadata_sample_names)]
        if any((':' in s) or (';' in s) for s in sample_ids):
            raise Unsupported("a sample name contains ':' or ';'")
        is_x = (args.chromosome == 'X')
        COHORT = Cohort.from_metadata(samples_df, sample_ids, is_x, args.x_par1_to, args.x_par2_from)
        rank = {name: r for r, name in enumerate(sample_ids)}
        record_parser = RecordParser(len(columns) - 9, [rank.get(name, -1) for name in columns[9:]])
        CHUNK_ROWS = pandas_chunk_rows(len(columns))
        SETTINGS = Settings(args.chromosome, args.force_rf, args.min_freq, args.min_hap_count, args.ignore_UTR,
                            args.skip_start_lost, len(args.output_cdna_fasta) > 0, selected_fixes(args.fix_bugs))

        print (('Chr ' + args.chromosome + ':'), "Reading", args.cdnas_fasta)
        cdnas = read_sequences(args.cdnas_fasta, set(ids))
        TRANSCRIPTS = [transcript_info(annotations_db, feature, cdnas) for feature in all_transcripts]
        del cdnas
        gc.collect()

        # the workers inherit everything above (and would also inherit unwritten output in stdout's buffer)
        sys.stdout.flush()
        pool = Pool(args.threads)
        state['pool'] = pool
        results = queue.Queue()
        if not fused:
            reader = VCFReader(input_vcf)
            state['reader'] = reader
            first_line_number, columns_again, first_line = read_header(reader)
            assert columns_again == columns
        outputs = Outputs(args, args.haplo_id_prefix)
        state['outputs'] = outputs

        print (('Chr ' + args.chromosome + ':'), 'Assigning variants to transcripts.')
        print (('Chr ' + args.chromosome + ':'), 'Computing haplotypes and protein haplotypes as the VCF is read, with',
               args.threads, 'worker(s).')
        done = {}
        in_flight = 0
        next_index = 0
        limit = IN_FLIGHT_PER_WORKER * max(args.threads, 1)

        def collect(block):
            nonlocal in_flight
            result = results.get(block=block)
            in_flight -= 1
            if isinstance(result, BaseException):
                raise result
            if result[0] == 'unsupported':
                raise Unsupported('%s: %s' % (TRANSCRIPTS[result[1]].id, result[2]))
            done[result[1]] = result[2:]

        def write_ready():
            nonlocal next_index
            while next_index in done:
                outputs.write(TRANSCRIPTS[next_index].id, *done.pop(next_index))
                next_index += 1

        for index, records in assign_records(reader, first_line, first_line_number, TRANSCRIPTS, record_parser):
            if records:
                while in_flight >= limit:
                    collect(True)
                pool.apply_async(haplotypes_of_transcript, (index, records), callback=results.put,
                                 error_callback=results.put)
                in_flight += 1
            else:
                if TRANSCRIPTS[index].cdna is None:
                    print('Transcript', TRANSCRIPTS[index].id, 'not in cDNA database, skipping!')
                done[index] = (1, 0, '', [], [], [])
            while not results.empty():
                collect(False)
            write_ready()
        reader.finish()
        state['reader'] = None
        while in_flight > 0:
            collect(True)
            write_ready()
        write_ready()
        assert next_index == len(TRANSCRIPTS) and not done

        pool.close()
        pool.join()
        state['pool'] = None
        print (('Chr ' + args.chromosome + ':'), 'Storing the result metadata:', args.output_file)
        print (('Chr ' + args.chromosome + ':'), 'Writing FASTA file:', args.output_fasta)
        if len(args.output_cdna_fasta) > 0:
            print (('Chr ' + args.chromosome + ':'), 'Writing cDNA FASTA file:', args.output_cdna_fasta)
        outputs.close()
        state['outputs'] = None
        print (('Chr ' + args.chromosome + ':'), "Done.")

    except Unsupported as e:
        if state['pool'] is not None:
            state['pool'].terminate()
            state['pool'].join()
        if state['reader'] is not None:
            state['reader'].close()
        if state['outputs'] is not None:
            state['outputs'].discard()
        if args.fix_bugs:
            # the reference has none of the fixes: its output would be that of -fix_bugs 0
            cannot_fall_back(__file__, str(e))
        argv = without_fix_bugs(sys.argv[1:])
        if fused:
            to_legacy_fused(__file__, args, str(e), argv)
        to_legacy(__file__, str(e), argv)
    except BaseException:
        if state['pool'] is not None:
            state['pool'].terminate()
        if state['reader'] is not None:
            state['reader'].close()
        raise


if __name__ == '__main__':
    parsed = parser.parse_args()
    if (parsed.input_vcf is None) == (parsed.filter_glob is None):
        parser.error('give either -i or -filter_glob')
    run(parsed)
