#!/usr/bin/env python3
"""Time fasta_simplify_headers.py (the workflow's last step) on simulated inputs with a given number of haplotype IDs.

    python3 bench/bench_simplify.py --ids 1000000,4000000 [--refs DIR] [--src src] [--work DIR]

The step reads the final FASTA, the haplotype table (without samples) and the whole-genome annotation DB, looks up the
transcript of every haplotype ID in the table (`haplo_df.loc[ID]`) and the gene of every transcript in the DB, and
writes the simplified FASTA and the header table. --src legacy/src measures the reference, which holds the whole
FASTA and table in memory; src/ streams them (IMPLEMENTATION_NOTES section 4.4).

For each --ids N this writes, into --work/N:
- final.fa: the MANE Select reference proteins (`generic_ensref|prot_<hex>|position_within_protein:0 start:0
  matching_proteins:ENST... reading_frame:-`, as merge_duplicate_seq.py and remove_UTR_seq.py leave them) and N / 1.5
  haplotype records (`generic_enshap|...` with one or two haplotype IDs each, alternately; the sequence of their
  transcript's protein), about the ratio of IDs to records of prohap.py's output at 500,000 samples (117k rows, 79k
  records per benchmark block);
- haplo.tsv.gz: the haplotype table as haplo_extract_sample_names.py writes it (20 columns), one row per ID.

Then it runs the script once, with PYTHONHASHSEED=0, and prints wall time, CPU time (rusage of the process from
wait4) and peak resident memory. Inputs are kept and reused.
"""
import argparse
import gzip
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COLUMNS = ['TranscriptID', 'chromosome', 'transcript_biotype', 'HaplotypeID', 'VCF_IDs', 'DNA_changes',
           'allele_frequencies', 'cDNA_changes', 'all_protein_changes', 'variant_types', 'protein_changes',
           'reading_frame', 'protein_prefix_length', 'start_missing', 'start_lost', 'splice_sites_affected',
           'occurrence_count', 'frequency', 'frequency_population', 'frequency_superpopulation']


def default_refs():
    root = os.environ.get('PROHAP_REFS', '')
    found = sorted(d for d in (os.listdir(root) if root and os.path.isdir(root) else []) if d.startswith('ensembl.v'))
    return os.path.join(root, found[0]) if len(found) == 1 else None


def read_fasta(fh):
    header, parts = None, []
    for line in fh:
        if line.startswith('>'):
            if header is not None:
                yield header, ''.join(parts)
            header, parts = line[1:].rstrip('\n'), []
        else:
            parts.append(line.strip())
    if header is not None:
        yield header, ''.join(parts)


def write_inputs(refs, n_ids, out):
    with open(os.path.join(refs, 'release.txt')) as fh:
        release = fh.read().strip()
    with open(os.path.join(refs, 'transcripts_reference_%s.csv' % release)) as fh:
        mane = {line.split(',')[1].strip() for line in list(fh)[1:]}
    reference = []          # (transcript, chromosome, sequence)
    with open(os.path.join(refs, 'ensembl_reference_proteinDB_%s_tagged.fa' % release)) as fh:
        for header, seq in read_fasta(fh):
            _, accession, description = header.split('|', 2)
            if accession in mane and '*' not in seq:
                chrom = description.split('chromosome:', 1)[1].split(':')[1] if 'chromosome:' in description else '1'
                reference.append((accession, chrom, seq))
    counters = {}
    record = 0
    ids = 0
    with open(os.path.join(out, 'final.fa'), 'w') as fasta, \
            gzip.open(os.path.join(out, 'haplo.tsv.gz'), 'wt', compresslevel=1) as table:
        table.write('\t'.join(COLUMNS) + '\n')
        for transcript, _, seq in reference:
            fasta.write('>generic_ensref|prot_%x|position_within_protein:0 start:0 matching_proteins:%s '
                        'reading_frame:-\n%s\n' % (record, transcript, seq))
            record += 1
        k = 0
        while ids < n_ids:
            transcript, chrom, seq = reference[k % len(reference)]
            haplotypes = []
            for _ in range(1 + k % 2):
                counters[chrom] = counters.get(chrom, 0) + 1
                hap = 'haplo_chr%s_%x' % (chrom, counters[chrom])
                haplotypes.append(hap)
                table.write('\t'.join([transcript, chrom, 'protein_coding', hap, 'rs%d' % ids, '100:A>G', '0.01',
                                       '10:A>G', '4:V>4:D', 'SAV', '4:V>4:D', '0', '0', 'False', 'False', '-', '3',
                                       '3e-06', 'GBR:0.001', 'EUR:0.001']) + '\n')
                ids += 1
            fasta.write('>generic_enshap|prot_%x|position_within_protein:0 start:0 matching_proteins:%s '
                        'reading_frame:%s\n%s\n' % (record, ','.join(haplotypes), ','.join('0' * len(haplotypes)),
                                                    seq))
            record += 1
            k += 1
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--ids', required=True, help='haplotype IDs, comma-separated sizes')
    parser.add_argument('--refs', default=default_refs(), help='references directory (default: from $PROHAP_REFS)')
    parser.add_argument('--src', default=os.path.join(ROOT, 'src'), help='script directory (default src/)')
    parser.add_argument('--work', default=os.path.join(ROOT, 'bench', 'data', 'simplify'),
                        help='input and output directory (default bench/data/simplify)')
    args = parser.parse_args(argv)
    if not args.refs:
        parser.error('--refs is required ($PROHAP_REFS does not hold exactly one ensembl.v<release>)')
    refs = os.path.abspath(os.path.expanduser(args.refs))
    args.src = os.path.abspath(args.src)        # the script runs in --work
    with open(os.path.join(refs, 'release.txt')) as fh:
        release = fh.read().strip()
    db = os.path.join(refs, 'Homo_sapiens.GRCh38.%s.chr_patch_hapl_scaff.db' % release)
    print('ids\trecords\twall_s\tcpu_s\tmax_mb')
    for n_ids in [int(x) for x in args.ids.split(',')]:
        out = os.path.join(args.work, str(n_ids))
        os.makedirs(out, exist_ok=True)
        done = os.path.join(out, 'records.txt')
        if not os.path.exists(done):
            records = write_inputs(refs, n_ids, out)
            with open(done, 'w') as fh:
                fh.write('%d\n' % records)
        with open(done) as fh:
            records = int(fh.read())
        command = [sys.executable, os.path.join(args.src, 'fasta_simplify_headers.py'), '-i', 'final.fa', '-hap_tsv',
                   'haplo.tsv.gz', '-db', db, '-o', 'final_simplified.fasta', '-header', 'final_header.tsv']
        start = time.time()
        proc = subprocess.Popen(command, cwd=out, env=dict(os.environ, PYTHONHASHSEED='0'))
        _, status, usage = os.wait4(proc.pid, 0)
        wall = time.time() - start
        if status != 0:
            sys.exit('fasta_simplify_headers.py failed (status %d) in %s' % (status, out))
        scale = 1 if sys.platform != 'darwin' else 1024     # ru_maxrss: KB on Linux, bytes on macOS
        print('%d\t%d\t%.1f\t%.1f\t%.0f' % (n_ids, records, wall, usage.ru_utime + usage.ru_stime,
                                           usage.ru_maxrss / scale / 1024), flush=True)


if __name__ == '__main__':
    main()
