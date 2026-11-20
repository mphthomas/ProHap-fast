from modules.common import read_fasta
from modules.fasta_stream import Accessions, records
import argparse
import gzip
import os
import pickle
import shutil
import tempfile

# ProHap-fast (IMPLEMENTATION_NOTES 1, D3, and 8.3, D4): the same output as the reference
# (legacy/src/merge_duplicate_seq.py), which keeps every record in memory, merges records with equal hash(sequence)
# and writes the merged entries in the order of their hash. Here the records are read one at a time and written,
# reduced to the fields the output needs, into PARTITIONS temporary files by the top bits of their hash, so that
# partition k holds the hashes of the k-th range in input order. Then each partition is read back, merged and
# sorted by hash in turn, and written with a running entry number: the same entries in the same order, holding one
# partition in memory instead of every record. If an accession repeats, read_fasta's dict would have replaced the
# earlier record; the run is then done again the dict way (merge_in_memory, D3), which is exact but holds everything.

PARTITIONS = 128            # at 25M records about 200,000 records each
BATCH = 500                 # records pickled at a time per partition (128 x 500 records buffered)

parser = argparse.ArgumentParser(description='Reads a FASTA file, aggregated all the duplicate sequences into one entry.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input FASTA file", metavar="FILE")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output FASTA file", metavar="FILE")

args = parser.parse_args()


def fields(protein):
	"""The reference's per-record fields: (matching_proteins, rfs, protein_start, seq_position)."""
	if ('hap' in protein['tag']) or ('var' in protein['tag']):
		matching_proteins = protein['description'].split('matching_proteins:', 1)[1].split(maxsplit=1)[0].split(';')
		rfs = protein['description'].split('reading_frame:', 1)[1].split(maxsplit=1)[0].split(';')
		protein_start = protein['description'].split('start:', 1)[1].split(maxsplit=1)[0]
		if ('position_within_protein' in protein['description']):
			seq_position = protein['description'].split('position_within_protein:', 1)[1].split(maxsplit=1)[0]
		else:
			seq_position = '0'
	elif 'ref' in protein['tag']:
		matching_proteins = protein['description'].split('matching_proteins:', 1)[1].split(maxsplit=1)[0].split(';')
		rfs = ['-']
		protein_start = '0'
		if ('position_within_protein' in protein['description']):
			seq_position = protein['description'].split('position_within_protein:', 1)[1].split(maxsplit=1)[0]
		else:
			seq_position = '0'
	else:
		matching_proteins = [protein['accession']]
		rfs = ['-']
		protein_start = '0'
		seq_position = '0'
	return matching_proteins, rfs, protein_start, seq_position


def add(result_by_hash, seq_hash, protein, matching_proteins, rfs, protein_start, seq_position):
	entry = result_by_hash.get(seq_hash)
	if entry is not None:
		entry['matching_proteins'].append(matching_proteins)
		entry['rfs'].append(rfs)
		entry['split_sequences'].append(protein['accession'])
		entry['tags'].append(protein['tag'])
		entry['start'].append(protein_start)
		entry['seq_position'].append(seq_position)
	else:
		result_by_hash[seq_hash] = {'hash': seq_hash, 'matching_proteins': [matching_proteins], 'split_sequences': [protein['accession']], 'tags': [protein['tag']], 'sequence': protein['sequence'], 'start': [protein_start], 'seq_position': [seq_position], 'rfs': [rfs]}


def write_entries(outfile, result_proteins, first_index):
	for i,protein in enumerate(result_proteins, first_index):
		matching_proteins = [ ','.join(plist) for plist in protein['matching_proteins'] ]
		rfs = [ ','.join(rflist) for rflist in protein['rfs'] ]
		description = 'position_within_protein:' + ';'.join(protein['seq_position']) + ' start:' + ';'.join(protein['start']) + ' matching_proteins:' + ';'.join(matching_proteins) + ' reading_frame:' + ';'.join(rfs) + ' split_sequences:' + ';'.join(protein['split_sequences'])

		tag = 'generic_'
		if ('generic_cont' in protein['tags']) or ('generic_sp' in protein['tags']):
			tag += 'cont'
		elif ('generic_ensref' in protein['tags']):
			tag += 'ensref'
		elif ('generic_ensvar' in protein['tags']):
			tag += 'ensvar'
		elif ('generic_var' in protein['tags']):
			tag += 'var'
		elif ('generic_manual' in protein['tags']):		# in case some sequences were added to the fasta manually
			tag += 'manual'
		elif ('generic_enshap' in protein['tags']):
			tag += 'enshap'
		elif ('generic_decoyvar' in protein['tags']):
			tag += 'decoyvar'
		else:
			tag += 'other'

		outfile.write('>' + tag + '|prot_' + hex(i)[2:] + '|' + description.replace('\n', '') + '\n')
		outfile.write(protein['sequence'] + '\n')


def open_output():
	return gzip.open(args.output_file, 'wt', compresslevel=1) if args.output_file.endswith('.gz') else open(args.output_file, 'w')


def merge_in_memory():
	"""D3: every record in a dict keyed by hash, sorted once."""
	all_proteins = read_fasta(args.input_file)
	result_by_hash = {}
	for protein in all_proteins.values():
		add(result_by_hash, hash(protein['sequence']), protein, *fields(protein))
	del all_proteins
	with open_output() as outfile:
		write_entries(outfile, [ result_by_hash[h] for h in sorted(result_by_hash) ], 0)


def partition_of(seq_hash):
	"""Partitions in the order of the hashes: the top bits of hash + 2**63."""
	return (seq_hash + (1 << 63)) * PARTITIONS >> 64


def merge_by_partitions(tmp_dir):
	"""Returns False, having written nothing, when an accession repeats."""
	paths = [os.path.join(tmp_dir, '%03d.pkl.gz' % k) for k in range(PARTITIONS)]
	files = [None] * PARTITIONS
	batches = [[] for _ in range(PARTITIONS)]
	accessions = Accessions()

	def flush(k):
		if files[k] is None:
			files[k] = gzip.open(paths[k], 'wb', compresslevel=1)
		pickle.dump(batches[k], files[k], protocol=pickle.HIGHEST_PROTOCOL)
		batches[k] = []

	for element_id, protein in records(args.input_file):
		accessions.add(element_id)
		seq_hash = hash(protein['sequence'])
		k = partition_of(seq_hash)
		batches[k].append((seq_hash, protein['tag'], protein['accession'], protein['sequence']) + fields(protein))
		if len(batches[k]) >= BATCH:
			flush(k)
	for k in range(PARTITIONS):
		if batches[k]:
			flush(k)
		if files[k] is not None:
			files[k].close()
	if accessions.repeated():
		return False

	written = 0
	with open_output() as outfile:
		for k in range(PARTITIONS):
			if files[k] is None:
				continue
			result_by_hash = {}
			with gzip.open(paths[k], 'rb') as fh:
				while True:
					try:
						batch = pickle.load(fh)
					except EOFError:
						break
					for seq_hash, tag, accession, sequence, matching_proteins, rfs, protein_start, seq_position in batch:
						protein = {'tag': tag, 'accession': accession, 'sequence': sequence}
						add(result_by_hash, seq_hash, protein, matching_proteins, rfs, protein_start, seq_position)
			os.remove(paths[k])
			entries = [ result_by_hash[h] for h in sorted(result_by_hash) ]
			del result_by_hash
			write_entries(outfile, entries, written)
			written += len(entries)
	return True


tmp_dir = tempfile.mkdtemp(prefix='merge_duplicate_seq.', dir=os.path.dirname(os.path.abspath(args.output_file)))
try:
	if not merge_by_partitions(tmp_dir):
		print('merge_duplicate_seq.py: an accession repeats in ' + args.input_file + '; merging in memory, as the reference reads it')
		merge_in_memory()
finally:
	shutil.rmtree(tmp_dir, ignore_errors=True)
