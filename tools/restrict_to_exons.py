"""Keep only the phased VCF records that overlap an exon of the listed transcripts (config `restrict_to_exons`).

    python3 tools/restrict_to_exons.py -i 'GLOB' -chr C -db CHR_DB -transcripts LIST -o OUT.vcf.gz [-threads N]

Input preparation, run the same way for both arms before the filter (IMPLEMENTATION_NOTES section 8.1, F2). The
reference reads every record of the chromosome twice (filter, then prohap.py) but uses only those whose REF overlaps
an exon of a listed transcript; with whole-genome input that is about 1 record in 50. This reads only the index
blocks that hold such records.

Every chunk of GLOB (in alphabetical order, as the filter reads them) must be indexed: bgzip-compressed VCF with a
.tbi or .csi index (`tabix -p vcf`), or BCF with a .csi index (`bcftools index`). VCF records are copied as they are
(`tabix -u -R`: a record that overlaps several regions is written once); BCF is converted to text by `bcftools view
-R --regions-overlap record`. OUT holds the first chunk's header and then the kept records of every chunk in order,
bgzip-compressed. A record is kept when [POS, POS + len(REF) - 1] meets an exon, the rule prohap.py uses, so every
record any transcript gets is kept, in the same order.

What changes: records outside the exons are gone, so the line numbers of the filtered VCF change, and with them the
IDs that prohap.py gives records whose ID is '.' (hex(line number)). Set IDs first (`bcftools annotate --set-id`) if
results must match a run without this option. Records of other contigs than the chromosome's are dropped too.
"""
import argparse
import glob
import os
import subprocess
import sys
import tempfile

import gffutils
import pandas as pd


class RestrictError(Exception):
    pass


def transcript_ids(path, chromosome):
    """The transcripts prohap.py takes from the list for this chromosome (before its start-codon check), as it reads
    them."""
    df = pd.read_csv(path)
    df['chromosome'] = df['chromosome'].apply(lambda x: str(x))
    return df[df['chromosome'] == chromosome]['transcriptID'].tolist()


def exon_regions(db_path, ids):
    """Sorted, merged [start, end] (1-based, inclusive) of the exons of the transcripts `ids`. A transcript missing
    from the database is skipped (prohap.py stops on it)."""
    db = gffutils.FeatureDB(db_path)
    exons = []
    for tid in ids:
        try:
            feature = db[tid]
        except gffutils.exceptions.FeatureNotFoundError:
            continue
        exons.extend((e.start, e.end) for e in db.children(feature, featuretype='exon'))
    merged = []
    for start, end in sorted(exons):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(r) for r in merged]


def contig_names(chromosome):
    names = {chromosome, 'chr' + chromosome}
    if chromosome.startswith('chr'):
        names.add(chromosome[3:])
    if chromosome in ('MT', 'chrMT', 'M', 'chrM'):
        names.update(('MT', 'chrMT', 'M', 'chrM'))
    return names


def index_of(path):
    for ext in ('.tbi', '.csi'):
        if os.path.isfile(path + ext):
            return path + ext
    return None


def listed_contigs(path, is_bcf):
    cmd = ['bcftools', 'index', '--stats', path] if is_bcf else ['tabix', '-l', path]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RestrictError('cannot read the index of %s: %s' % (path, proc.stderr.strip()))
    return [line.split('\t')[0] for line in proc.stdout.splitlines() if line]


def chunk_command(path, chromosome, regions, first):
    """The command that writes the kept records of one chunk as VCF text (with the header when `first`); None when
    the chunk has no record of the chromosome and is not the first."""
    is_bcf = path.endswith('.bcf')
    if index_of(path) is None:
        raise RestrictError(
            '%s has no index (.tbi or .csi). restrict_to_exons reads only the indexed blocks that hold exons, so every '
            'VCF chunk must be bgzip-compressed and indexed (bgzip FILE; tabix -p vcf FILE.gz), and BCF indexed '
            '(bcftools index FILE). Index the chunks or set restrict_to_exons: False.' % path)
    contigs = listed_contigs(path, is_bcf)
    wanted = [c for c in contigs if c in contig_names(chromosome)]
    if len(wanted) > 1:
        raise RestrictError('%s has records under more than one name of chromosome %s: %s'
                            % (path, chromosome, ', '.join(wanted)))
    others = [c for c in contigs if c not in wanted]
    if others:
        print('restrict_to_exons: %s also holds contig(s) %s; their records are dropped'
              % (path, ', '.join(others[:5]) + (', ...' if len(others) > 5 else '')), file=sys.stderr)
    if not wanted or not regions:
        if not first:
            return None
        return ['bcftools', 'view', '--no-version', '-h', '-Ov', path] if is_bcf else ['tabix', '-H', path]
    contig = wanted[0]
    if is_bcf:
        return (['bcftools', 'view', '--no-version', '-Ov', '--regions-overlap', 'record', '-R', regions_file(contig, regions)]
                + ([] if first else ['-H']) + [path])
    return ['tabix'] + (['-h'] if first else []) + ['-u', '-R', regions_file(contig, regions), path]


_regions_files = {}
_tmp = None


def regions_file(contig, regions):
    """A tab-separated CHROM, START, END file (1-based, inclusive) of `regions` on `contig`."""
    if contig not in _regions_files:
        path = os.path.join(_tmp.name, 'regions_%d.tsv' % len(_regions_files))
        with open(path, 'w') as fh:
            fh.write(''.join('%s\t%d\t%d\n' % (contig, s, e) for s, e in regions))
        _regions_files[contig] = path
    return _regions_files[contig]


def restrict(files, chromosome, regions, output, threads=1):
    """Write the records of `files` (in this order) that overlap `regions` to `output` (bgzip)."""
    global _tmp
    if not files:
        raise RestrictError('no VCF matches the input pattern for chromosome %s' % chromosome)
    _regions_files.clear()
    _tmp = tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(output)))
    try:
        commands = [chunk_command(path, chromosome, regions, i == 0) for i, path in enumerate(files)]
        with open(output + '.tmp', 'wb') as out:
            bgzip = subprocess.Popen(['bgzip', '-c', '-l', '1', '-@', str(max(1, threads))],
                                     stdin=subprocess.PIPE, stdout=out)
            try:
                for path, cmd in zip(files, commands):
                    if cmd is None:
                        continue
                    proc = subprocess.run(cmd, stdout=bgzip.stdin, stderr=subprocess.PIPE)
                    if proc.returncode != 0:
                        raise RestrictError('%s failed on %s: %s'
                                            % (cmd[0], path, proc.stderr.decode(errors='replace').strip()))
            finally:
                bgzip.stdin.close()
                code = bgzip.wait()
            if code != 0:
                raise RestrictError('bgzip failed writing ' + output)
        os.replace(output + '.tmp', output)
    finally:
        _tmp.cleanup()
        if os.path.exists(output + '.tmp'):
            os.remove(output + '.tmp')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('-i', dest='input_glob', required=True, help='the VCF chunks of the chromosome (glob)')
    p.add_argument('-chr', dest='chromosome', required=True)
    p.add_argument('-db', dest='db', required=True, help='gffutils database of the chromosome')
    p.add_argument('-transcripts', dest='transcripts', required=True, help='the transcript list (CSV)')
    p.add_argument('-o', dest='output', required=True, help='output VCF (bgzip)')
    p.add_argument('-threads', dest='threads', type=int, default=1)
    args = p.parse_args(argv)
    # index files cannot be chunks (a pattern such as chr22* matches them)
    files = [f for f in sorted(glob.glob(args.input_glob)) if not f.endswith(('.tbi', '.csi'))]
    regions = exon_regions(args.db, transcript_ids(args.transcripts, args.chromosome))
    try:
        restrict(files, args.chromosome, regions, args.output, args.threads)
    except RestrictError as e:
        print('restrict_to_exons: error: %s' % e, file=sys.stderr)
        return 1
    print('Chr %s: kept the records on %d exon regions (%d bp) of %d VCF chunk(s)'
          % (args.chromosome, len(regions), sum(e - s + 1 for s, e in regions), len(files)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
