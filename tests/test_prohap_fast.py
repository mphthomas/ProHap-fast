"""src/prohap.py (R1c) against legacy/src/prohap.py, as scripts, on a crafted chromosome (a gffutils DB built from a
GTF, a cDNA FASTA, sample metadata and a VCF; see fast_cases.make_protein_case): identical outputs with the fast path
taken, and identical outputs (or both runs stopping) for every input that hands the run to the reference.

    pytest tests/test_prohap_fast.py
"""
import concurrent.futures
import gzip
import os
import subprocess
import sys
from pathlib import Path

import pytest

gffutils = pytest.importorskip('gffutils')

from fast_cases import make_protein_case  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
TREES = {'legacy': ROOT / 'legacy' / 'src', 'fast': ROOT / 'src'}
FALLBACK = 'ProHap-fast: falling back to the reference code'
OUTPUTS = ('haplo.tsv.gz', 'haplo.fa', 'cdna.fa', 'prohap.log')


def write_case(path, seed=1, vcf_edit=None, metadata_edit=None, list_edit=None, gz=True):
    features, db, cdnas, samples, text = make_protein_case(seed, n_samples=30)
    gtf = []
    for f in features:
        attrs = 'gene_id "G%s"; transcript_id "%s";' % (f.id, f.id)
        if f.biotype is not None:
            attrs += ' transcript_biotype "%s";' % f.biotype
        gtf.append('\t'.join(['22', 'test', 'transcript', str(f.start), str(f.end), '.', f.strand, '.', attrs]))
        for kind, spans in sorted(db.data[f.id].items()):
            for s, e in spans:
                gtf.append('\t'.join(['22', 'test', kind, str(s), str(e), '.', f.strand, '.', attrs]))
    (path / 'in.gtf').write_text('\n'.join(gtf) + '\n')
    gffutils.create_db(str(path / 'in.gtf'), dbfn=str(path / 'in.db'), force=True, keep_order=True,
                       merge_strategy='merge', sort_attribute_values=True, disable_infer_genes=True,
                       disable_infer_transcripts=True)
    with open(path / 'cdna.fa', 'w') as out:
        for tid, entry in cdnas.items():
            seq = entry['sequence']
            out.write('>%s.1 cdna chromosome:test\n' % tid + ''.join(seq[i:i + 60] + '\n'
                                                                    for i in range(0, len(seq), 60)))
    listed = 'chromosome,transcriptID\n' + ''.join('22,%s\n' % f.id for f in features) + '21,ENST_OTHER\n'
    if list_edit:
        listed = list_edit(listed)
    (path / 'transcripts.csv').write_text(listed)
    metadata = 'Sample name\tSex\tPopulation code\tSuperpopulation code\n' + ''.join(
        '%s\t%s\t%s\t%s\n' % (s, ['male', 'female'][i % 2], ['GBR', 'YRI', 'CHB'][i % 3], ['EUR', 'AFR', 'EAS'][i % 3])
        for i, s in enumerate(samples[:-2])) + 'M01\tmale\tGBR\tEUR\n'
    if metadata_edit:
        metadata = metadata_edit(metadata)
    (path / 'samples.tsv').write_text(metadata)
    if vcf_edit:
        text = vcf_edit(text)
    data = text.encode('utf-8') if isinstance(text, str) else text
    vcf = path / ('in.vcf.gz' if gz else 'in.vcf')
    if gz:
        with gzip.open(vcf, 'wb') as f:
            f.write(data)
    else:
        vcf.write_bytes(data)
    return vcf


def run(arm, case, vcf, extra=(), threads=2, chromosome='22'):
    out = case / arm
    out.mkdir()
    (out / 'tmp').mkdir()
    command = [sys.executable, str(ROOT / 'tools' / 'run_fork.py'), str(TREES[arm] / 'prohap.py'), '-i', str(vcf),
               '-db', str(case / 'in.db'), '-transcripts', str(case / 'transcripts.csv'), '-cdna',
               str(case / 'cdna.fa'), '-s', str(case / 'samples.tsv'), '-chr', chromosome, '-min_hap_freq', '-1',
               '-min_hap_count', '0', '-acc_prefix', 'enshap_22', '-id_prefix', 'haplo_chr22', '-threads',
               str(threads), '-log', 'prohap.log', '-tmp_dir', 'tmp', '-output_csv', 'haplo.tsv.gz', '-output_fasta',
               'haplo.fa', '-output_cdna_fasta', 'cdna.fa'] + list(extra)
    proc = subprocess.run(command, cwd=out, capture_output=True, text=True,
                          env=dict(os.environ, PYTHONHASHSEED='0', PYTHONNOUSERSITE='1'))
    files = {}
    for name in OUTPUTS:
        if (out / name).exists():
            data = (out / name).read_bytes()
            files[name] = gzip.decompress(data) if data[:2] == b'\x1f\x8b' else data
    return proc, files


def run_both(case, vcf, extra=(), threads=2, chromosome='22'):
    legacy = run('legacy', case, vcf, extra, threads, chromosome)
    fast = run('fast', case, vcf, extra, threads, chromosome)
    return legacy, fast


def first_genotype_line(text):
    lines = text.split('\n')
    i = next(k for k, line in enumerate(lines) if line and not line.startswith('#'))
    return lines, i


def edit_token(old, new):
    def edit(text):
        lines, i = first_genotype_line(text)
        fields = lines[i].split('\t')
        k = next(j for j in range(9, len(fields)) if fields[j] == old)
        fields[k] = new
        lines[i] = '\t'.join(fields)
        return '\n'.join(lines)
    return edit


def edit_field(index, value):
    def edit(text):
        lines, i = first_genotype_line(text)
        fields = lines[i].split('\t')
        fields[index] = value
        lines[i] = '\t'.join(fields)
        return '\n'.join(lines)
    return edit


def duplicate_sample(text):
    head, rest = text.split('\n#CHROM')
    header, body = rest.split('\n', 1)
    names = header.split('\t')
    names[-1] = names[-2]
    return head + '\n#CHROM' + '\t'.join(names) + '\n' + body


def empty_vcf(text):
    head, rest = text.split('\n#CHROM')
    return head + '\n#CHROM' + rest.split('\n')[0] + '\n'


NO_UTR = ('-ignore_UTR', '0')       # the crafted reading frames meet a stop codon early; keep the UTR variants
# name: (write_case arguments, run_both arguments, fast path?, does the reference stop?)
CASES = {
    # the fast path
    'threads1': (dict(seed=1), dict(threads=1), True, False),
    'threads3_plain_vcf': (dict(seed=2, gz=False), dict(threads=3, extra=NO_UTR), True, False),
    'threads2': (dict(seed=3), dict(extra=NO_UTR), True, False),
    'min_hap_freq': (dict(seed=4), dict(extra=('-min_hap_freq', '0.05')), True, False),
    'min_hap_count': (dict(seed=4), dict(extra=('-min_hap_count', '3') + NO_UTR), True, False),
    'start_lost_and_no_start': (dict(seed=4), dict(extra=('-skip_start_lost', '0', '-require_start', '0')), True,
                                False),
    'empty_vcf': (dict(seed=1, vcf_edit=empty_vcf), {}, True, False),
    # handed to the reference, which writes the same output
    'unphased': (dict(seed=2, vcf_edit=edit_token('0|0', '0/1')), dict(extra=NO_UTR), False, False),
    'multi_character_allele': (dict(seed=2, vcf_edit=edit_token('0|0', '10|0')), dict(extra=NO_UTR), False, False),
    'extra_format_field': (dict(seed=2, vcf_edit=lambda t: t.replace('\tGT\t', '\tGT:DS\t')), dict(extra=NO_UTR),
                           False, False),
    'quote': (dict(seed=2, vcf_edit=edit_field(2, '"rs1"')), dict(extra=NO_UTR), False, False),
    'carriage_return': (dict(seed=2, vcf_edit=lambda t: t.replace('\n', '\r\n')), dict(extra=NO_UTR), False, False),
    'non_ascii': (dict(seed=2, vcf_edit=edit_field(7, 'AF=0.1;NOTE=é')), dict(extra=NO_UTR), False, False),
    'duplicate_sample': (dict(seed=2, vcf_edit=duplicate_sample), dict(extra=NO_UTR), False, False),
    'info_not_text': (dict(seed=2, vcf_edit=edit_field(7, '1')), dict(extra=NO_UTR), False, False),
    'population_with_colon': (dict(seed=2, metadata_edit=lambda m: m.replace('\tGBR\t', '\tGB:R\t')),
                              dict(extra=NO_UTR), False, False),
    'af_threshold': (dict(seed=2), dict(extra=('-af', '0.2') + NO_UTR), False, False),
    # handed to the reference, which stops
    'duplicate_transcript': (dict(seed=2, list_edit=lambda t: t + t.split('\n')[1] + '\n'), {}, False, True),
    'duplicate_metadata_row': (dict(seed=2, metadata_edit=lambda m: m + m.split('\n')[3] + '\n'), {}, False, True),
    'no_sample_in_metadata': (dict(seed=2, metadata_edit=lambda m: m.replace('S0', 'Z0')), {}, False, True),
    'no_transcript_for_the_records': (dict(seed=2, list_edit=lambda t: t.replace('21,ENST_OTHER\n', '')),
                                      dict(chromosome='21'), False, True),
}


@pytest.fixture(scope='module')
def results(tmp_path_factory):
    """Every case, run by both scripts, eight cases at a time."""
    root = tmp_path_factory.mktemp('cases')

    def one(name):
        case = root / name
        case.mkdir()
        write_args, run_args, _, _ = CASES[name]
        return name, run_both(case, write_case(case, **write_args), **run_args)
    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        return dict(pool.map(one, CASES))


@pytest.mark.parametrize('name', list(CASES))
def test_case(results, name):
    _, _, fast_path, stops = CASES[name]
    (legacy, legacy_files), (fast, fast_files) = results[name]
    assert (FALLBACK in fast.stderr) == (not fast_path), fast.stderr
    if stops:
        assert legacy.returncode != 0 and fast.returncode != 0, (legacy.stderr, fast.stderr)
        return
    assert legacy.returncode == 0, legacy.stderr
    assert fast.returncode == 0, fast.stderr
    assert fast_files == legacy_files
    if name == 'empty_vcf':
        assert set(fast_files) == {'haplo.tsv.gz', 'haplo.fa'}
    elif not name.startswith('min_hap'):
        assert fast_files['haplo.tsv.gz'].count(b'\n') > 10 and fast_files['prohap.log'].count(b'\n') > 3
