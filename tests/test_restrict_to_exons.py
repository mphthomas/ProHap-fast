"""tools/restrict_to_exons.py (F2, config restrict_to_exons).

The script must keep every record whose [POS, POS + len(REF) - 1] meets an exon of a listed transcript, in file
order, each once, after the first chunk's header: then prohap.py gives every transcript the same records as without
the option. Checked against a direct implementation of that rule on crafted chunks (long deletions across exons and
gaps, exon edges, several records at one position, records of another contig), and on the 1kG sample data: with
IDs set, a workflow run with the option writes the same outputs as one without (the filtered VCF aside).

    pytest tests/test_restrict_to_exons.py           # the workflow test needs PROHAP_REFS
"""
import gzip
import os
import random
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'tools'))
import restrict_to_exons as rte  # noqa: E402

pytestmark = pytest.mark.skipif(shutil.which('tabix') is None or shutil.which('bgzip') is None,
                                reason='needs tabix and bgzip (htslib)')

HEADER = ('##fileformat=VCFv4.2\n##contig=<ID=22>\n##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
          '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n')


def record(contig, pos, ref, alt='A', rid='.', info='AF=0.1'):
    return '%s\t%d\t%s\t%s\t%s\t.\tPASS\t%s\tGT\t0|1\t1|0\n' % (contig, pos, rid, ref, alt, info)


def write_chunk(path, records, header=HEADER):
    plain = str(path)[:-3]
    with open(plain, 'w') as fh:
        fh.write(header + ''.join(records))
    subprocess.run(['bgzip', '-f', plain], check=True)
    subprocess.run(['tabix', '-f', '-p', 'vcf', str(path)], check=True)
    return str(path)


def expected(chunks, regions, contigs=('22', 'chr22')):
    """The rule, directly: the first chunk's header, then every record that overlaps a region, in file order."""
    out = []
    for i, recs in enumerate(chunks):
        if i == 0:
            out.append(HEADER)
        for r in recs:
            f = r.split('\t')
            pos, ref = int(f[1]), f[3]
            if f[0] in contigs and any(pos <= e and pos + len(ref) - 1 >= s for s, e in regions):
                out.append(r)
    return ''.join(out)


def read(path):
    with gzip.open(path, 'rt') as fh:
        return fh.read()


def random_records(rng, start, n, contig='22'):
    pos = start
    recs = []
    for _ in range(n):
        pos += rng.choice([0, 0, 1, 3, 10, 40])
        ref = rng.choice(['A', 'C', 'GT', 'ACGTACGTAC', 'A' * rng.randint(2, 300)])
        recs.append(record(contig, pos, ref, rid='r%d' % len(recs)))
    return recs


def test_crafted_edges(tmp_path):
    regions = [(100, 200), (300, 300), (1000, 1100)]
    recs = [
        record('22', 50, 'A'),                  # before
        record('22', 60, 'A' * 41),             # ends at 100: kept
        record('22', 60, 'A' * 40),             # ends at 99: not
        record('22', 99, 'A'),
        record('22', 100, 'A'),
        record('22', 150, 'AC' * 100),          # 150-349: across the gap into 300-300, written once
        record('22', 150, 'A'),                 # same position
        record('22', 150, 'C', alt='G'),
        record('22', 200, 'A'),
        record('22', 201, 'A'),
        record('22', 250, 'A' * 900),           # 250-1149: meets all but the first, written once
        record('22', 299, 'A'),
        record('22', 300, 'A'),
        record('22', 301, 'A'),
        record('22', 1100, 'A'),
        record('22', 1101, 'A'),
    ]
    path = write_chunk(tmp_path / 'c.vcf.gz', recs)
    out = str(tmp_path / 'out.vcf.gz')
    rte.restrict([path], '22', regions, out)
    assert read(out) == expected([recs], regions)
    assert read(out).count(record('22', 250, 'A' * 900)) == 1


@pytest.mark.parametrize('seed', range(6))
def test_random_chunks(tmp_path, seed):
    rng = random.Random(seed)
    regions = []
    pos = rng.randint(1, 500)
    for _ in range(rng.randint(1, 40)):
        start = pos + rng.randint(2, 400)
        end = start + rng.randint(0, 200)
        regions.append((start, end))
        pos = end
    chunks, paths, start = [], [], 1
    for k in range(rng.randint(1, 3)):
        recs = random_records(rng, start, rng.randint(0, 400))
        start = (int(recs[-1].split('\t')[1]) + 1) if recs else start
        chunks.append(recs)
        paths.append(write_chunk(tmp_path / ('chr22.%d.vcf.gz' % k), recs))
    out = str(tmp_path / 'out.vcf.gz')
    rte.restrict(paths, '22', regions, out)
    assert read(out) == expected(chunks, regions)


def test_chr_prefix_and_other_contigs(tmp_path, capsys):
    regions = [(100, 200)]
    recs = [record('chr21', 150, 'A'), record('chr22', 120, 'A'), record('chr22', 500, 'A'), record('chrX', 150, 'A')]
    path = write_chunk(tmp_path / 'c.vcf.gz', recs)
    out = str(tmp_path / 'out.vcf.gz')
    rte.restrict([path], '22', regions, out)
    assert read(out) == HEADER + record('chr22', 120, 'A')
    assert 'chr21' in capsys.readouterr().err


def test_first_chunk_without_the_chromosome(tmp_path):
    regions = [(100, 200)]
    a = write_chunk(tmp_path / 'a.vcf.gz', [record('21', 150, 'A')])
    b = write_chunk(tmp_path / 'b.vcf.gz', [record('22', 150, 'A')], header=HEADER.replace('S1', 'OTHER'))
    out = str(tmp_path / 'out.vcf.gz')
    rte.restrict([a, b], '22', regions, out)
    assert read(out) == HEADER + record('22', 150, 'A')      # the first chunk's header, the second's records


def test_no_regions_keeps_the_header(tmp_path):
    path = write_chunk(tmp_path / 'c.vcf.gz', [record('22', 150, 'A')])
    out = str(tmp_path / 'out.vcf.gz')
    rte.restrict([path], '22', [], out)
    assert read(out) == HEADER


def test_unindexed_input_is_an_error(tmp_path):
    path = write_chunk(tmp_path / 'c.vcf.gz', [record('22', 150, 'A')])
    os.remove(path + '.tbi')
    with pytest.raises(rte.RestrictError, match='no index'):
        rte.restrict([path], '22', [(100, 200)], str(tmp_path / 'out.vcf.gz'))
    assert not os.path.exists(tmp_path / 'out.vcf.gz')


@pytest.mark.skipif(shutil.which('bcftools') is None, reason='needs bcftools')
def test_bcf(tmp_path):
    rng = random.Random(7)
    regions = [(200, 400), (900, 950), (2000, 2600)]
    recs = random_records(rng, 1, 300)
    vcf = write_chunk(tmp_path / 'c.vcf.gz', recs)
    bcf = str(tmp_path / 'c.bcf')
    subprocess.run(['bcftools', 'view', '--no-version', '-Ob', '-o', bcf, vcf], check=True)
    subprocess.run(['bcftools', 'index', bcf], check=True)
    out = str(tmp_path / 'out.vcf.gz')
    rte.restrict([bcf], '22', regions, out)
    kept = [line.split('\t')[2] for line in read(out).splitlines() if not line.startswith('#')]
    want = [line.split('\t')[2] for line in expected([recs], regions).splitlines() if not line.startswith('#')]
    assert kept == want


def test_exon_regions(tmp_path):
    gffutils = pytest.importorskip('gffutils')
    gtf = []
    for tid, exons in (('T1', [(100, 200), (500, 600)]), ('T2', [(150, 250), (601, 700)]), ('T3', [(900, 950)])):
        attrs = 'gene_id "G%s"; transcript_id "%s";' % (tid, tid)
        gtf.append('\t'.join(['22', 't', 'transcript', str(exons[0][0]), str(exons[-1][1]), '.', '+', '.', attrs]))
        gtf += ['\t'.join(['22', 't', 'exon', str(s), str(e), '.', '+', '.', attrs]) for s, e in exons]
    (tmp_path / 'a.gtf').write_text('\n'.join(gtf) + '\n')
    gffutils.create_db(str(tmp_path / 'a.gtf'), dbfn=str(tmp_path / 'a.db'), disable_infer_genes=True,
                       disable_infer_transcripts=True, merge_strategy='create_unique')
    (tmp_path / 'list.csv').write_text('chromosome,transcriptID\n22,T1\n22,T2\n21,T3\n22,MISSING\n')
    ids = rte.transcript_ids(str(tmp_path / 'list.csv'), '22')
    assert ids == ['T1', 'T2', 'MISSING']
    assert rte.exon_regions(str(tmp_path / 'a.db'), ids) == [(100, 250), (500, 700)]


# ------------------------------------ the workflow, with and without the option ------------------------------------

REFS = Path(os.environ.get('PROHAP_REFS') or '/opt/prohap/refs').expanduser()


@pytest.mark.skipif(not (REFS / 'ensembl.v113').is_dir(), reason='needs the release 113 references (PROHAP_REFS)')
@pytest.mark.skipif(shutil.which('snakemake') is None, reason='needs snakemake')
def test_workflow_same_outputs_with_ids_set(tmp_path):
    """1kG chromosomes 1 (two chunks), 21 and 22 with every ID set: the fast arm writes the same outputs with
    restrict_to_exons as without, and the legacy arm with it the same as the fast arm with it."""
    data = tmp_path / 'input'
    data.mkdir()
    names = ['sample_chr1_1kGP_GRCh38_phased.vcf', 'sample_chr1_1kGP_GRCh38_phased_2.vcf',
             'sample_chr21_1kGP_GRCh38_phased.vcf', 'sample_chr22_1kGP_GRCh38_phased.vcf']
    with tarfile.open(ROOT / 'sample_data' / 'sample_1kGP_common_global.tar.gz') as tar:
        for i, name in enumerate(names):
            text = tar.extractfile(name).read().decode()
            lines = text.splitlines(keepends=True)
            out = []
            for n, line in enumerate(lines):
                if not line.startswith('#'):
                    f = line.split('\t', 3)
                    line = '\t'.join(f[:2] + ['id%d_%d' % (i, n)] + f[3:])
                out.append(line)
            write_chunk(data / (name + '.gz'), out, header='')
    fixture = (ROOT / 'tests' / 'fixtures' / '1kg' / 'threshold0.yaml').read_text()
    fixture = fixture.replace('chromosomes: ["1", "21", "22", "X"]', 'chromosomes: ["1", "21", "22"]')
    fixture = fixture.replace('sample_chr{chr}_1kGP_GRCh38_phased*.vcf"', 'sample_chr{chr}_1kGP_GRCh38_phased*.vcf.gz"')
    assert 'phased*.vcf.gz"' in fixture and '"1", "21", "22"]' in fixture
    runs = {}
    for name, arm, restrict in (('plain', 'fast', False), ('fast', 'fast', True), ('legacy', 'legacy', True)):
        d = tmp_path / name
        d.mkdir()
        (d / 'config.yaml').write_text(fixture + '\nrestrict_to_exons: %s\n' % restrict)
        os.symlink(data, d / 'sample_data')
        for f in ('igsr_samples.tsv', 'crap.fasta'):
            os.symlink(ROOT / f, d / f)
        os.symlink(ROOT / 'tests' / 'fixtures' / '1kg' / 'transcripts.csv', d / 'transcripts.csv')
        env = dict(os.environ, PROHAP_REFS=str(REFS), PYTHONHASHSEED='0')
        proc = subprocess.run(['snakemake', '--snakefile', str(ROOT / 'Snakefile'), '--directory', str(d),
                               '--cores', '4', '--notemp', '--config', 'implementation=' + arm],
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
        assert proc.returncode == 0, proc.stdout[-3000:]
        runs[name] = d

    def outputs(d, skip):
        files = {}
        for p in sorted(d.rglob('*')):
            rel = str(p.relative_to(d))
            if p.is_file() and not p.is_symlink() and not rel.startswith('.snakemake') and rel != 'config.yaml' \
                    and not any(s in rel for s in skip):
                with open(p, 'rb') as fh:
                    raw = fh.read()
                files[rel] = gzip.decompress(raw) if raw[:2] == b'\x1f\x8b' else raw
        return files

    plain = outputs(runs['plain'], ['_phased_filtered.vcf.gz'])
    fast = outputs(runs['fast'], ['_phased_filtered.vcf.gz', '_exons.vcf.gz'])
    assert sorted(plain) == sorted(fast)
    assert [f for f in plain if plain[f] != fast[f]] == []
    assert len(plain['results/1kg_haplotypes.tsv'].splitlines()) > 100
    # fewer records in the restricted filtered VCF, and the arms agree on everything
    for c in ('1', '21', '22'):
        name = 'data/vcf/phased/chr%s_phased_filtered.vcf.gz' % c
        assert len(read(runs['fast'] / name)) < len(read(runs['plain'] / name))
    legacy = outputs(runs['legacy'], [])
    fast_all = outputs(runs['fast'], [])
    assert sorted(legacy) == sorted(fast_all)
    assert [f for f in legacy if legacy[f] != fast_all[f]] == []
