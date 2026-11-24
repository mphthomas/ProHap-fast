"""src/vcf_filter_fix.py (F1) against legacy/src/vcf_filter_fix.py, as scripts, on small crafted VCFs: one per case the
fast path handles and one per case it hands over (a line to the reference's own code, or the whole file to the
reference script):

    pytest tests/test_filter_fast.py
"""
import gzip
import os
import random
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
TREES = {'legacy': ROOT / 'legacy' / 'src', 'fast': ROOT / 'src'}
FALLBACK = 'ProHap-fast: falling back to the reference code'
HEADER = ['##fileformat=VCFv4.2', '##contig=<ID=22>']
COLUMNS = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']


def vcf(samples, records, header=HEADER):
    """records: (pos, id, ref, alt, info, format, genotypes) with genotypes a list of tokens."""
    lines = header + ['\t'.join(COLUMNS + samples)]
    for pos, vcf_id, ref, alt, info, fmt, genotypes in records:
        lines.append('\t'.join([str(22), str(pos), vcf_id, ref, alt, '.', 'PASS', info, fmt] + genotypes))
    return '\n'.join(lines) + '\n'


def write(path, text):
    data = text.encode('utf-8') if isinstance(text, str) else text
    if str(path).endswith('.gz'):
        with gzip.open(path, 'wb') as f:
            f.write(data)
    else:
        path.write_bytes(data)


def run_both(tmp_path, pattern, af='0', extra=()):
    results = {}
    for arm, src in TREES.items():
        out = tmp_path / ('out_%s.vcf.gz' % arm)
        proc = subprocess.run([sys.executable, str(src / 'vcf_filter_fix.py'), '-i', str(pattern), '-chr', '22',
                               '-af', af, '-o', str(out)] + list(extra), capture_output=True, text=True,
                              env=dict(os.environ, PYTHONHASHSEED='0'))
        output = gzip.decompress(out.read_bytes()) if out.exists() and proc.returncode == 0 else None
        results[arm] = proc.returncode, proc.stdout, output, proc.stderr
    return results


def check(results, fallback=False):
    legacy, fast = results['legacy'], results['fast']
    assert legacy[0] == 0, legacy[3]
    assert fast[:3] == legacy[:3]
    assert (FALLBACK in fast[3]) == fallback, fast[3]
    return legacy[2].decode()


def random_genotypes(rng, n, tokens=('0|0', '0|0', '0|0', '0|1', '1|0', '1|1')):
    return [rng.choice(tokens) for _ in range(n)]


SAMPLES = ['S%d' % i for i in range(1, 11)]


@pytest.mark.parametrize('af', ['0', '0.05', '-1'])
def test_biallelic_and_multiallelic(tmp_path, af):
    rng = random.Random(1)
    records = []
    for i in range(40):
        alt = rng.choice(['G', 'GT', 'C,G', 'A,C,T', 'C,G,T,AA,AC,AG,AT,CA,CC'])     # up to 9 ALT alleles
        n_alt = alt.count(',') + 1
        tokens = ['%d|%d' % (rng.randint(0, n_alt), rng.randint(0, n_alt)) for _ in SAMPLES]
        afs = ','.join(str(round(rng.random() / 5, 3)) for _ in range(n_alt))
        info = rng.choice(['AC=1;AN=20;AF=' + afs, 'AF=' + afs + ';AC=1', 'AC=1;AN=20;AF=' + afs + ';DP=7'])
        records.append((100 + i * 10, rng.choice(['.', 'rs%d' % i]), 'A', alt, info, 'GT', tokens))
    write(tmp_path / 'in.vcf', vcf(SAMPLES, records))
    out = check(run_both(tmp_path, tmp_path / 'in.vcf', af))
    assert 'MAF=' in out


def test_haploid_missing_and_unphased_tokens(tmp_path):
    """chrX males written haploid, missing alleles and unphased genotypes: the fast path handles them all."""
    rng = random.Random(2)
    records = []
    for i in range(30):
        tokens = random_genotypes(rng, 10, ('0', '1', '.', '0|0', '1|0', '0|1', '.|.', '0/1', '1/1', './0', '2|1'))
        records.append((100 + i, '.', 'A', rng.choice(['G', 'C,G', 'C,G,T']), 'AF=0.1,0.2,0.3', 'GT', tokens))
    records.append((500, '.', 'A', 'G', 'AF=0.1', 'GT', ['0'] * 10))                     # all haploid
    write(tmp_path / 'in.vcf.gz', vcf(SAMPLES, records))
    out = check(run_both(tmp_path, tmp_path / 'in.vcf.gz'))
    assert '1|0\t' in out and '.|0' in out and '0/1' in out


def test_chunks_headers_and_empty_files(tmp_path):
    """Chunks in alphabetical order, only the first file's header, a chunk with a header only, an empty file, a last
    line without a newline."""
    rng = random.Random(3)
    def records(start):
        return [(start + i, '.', 'C', 'T', 'AF=0.2', 'GT', random_genotypes(rng, 10)) for i in range(5)]
    write(tmp_path / 'c.1.vcf.gz', vcf(SAMPLES, records(100)))
    write(tmp_path / 'c.2.vcf.gz', vcf(SAMPLES, [], header=HEADER + ['##second=1']))
    write(tmp_path / 'c.3.vcf.gz', vcf(SAMPLES, records(200)).rstrip('\n'))
    write(tmp_path / 'c.4.vcf.gz', '')
    results = run_both(tmp_path, tmp_path / 'c.*.vcf.gz')
    out = check(results)
    assert out.count('#CHROM') == 1 and len(out.splitlines()) == 3 + 10
    assert results['fast'][1].count('is empty!') == 2


def test_no_input_files(tmp_path):
    check(run_both(tmp_path, tmp_path / 'none*.vcf'))


@pytest.mark.parametrize('af', ['0', '0.05'])
def test_format_fields_after_gt(tmp_path, af):
    """FORMAT with GT first and more fields (GT:DS, GT:AD:DP:GQ, GT:): the reference keeps the GT subfield of each
    token; so does the fast path. Includes trailing fields left out, whole-token '.', haploid tokens, unphased and
    missing genotypes, multi-allelic lines, and tokens the fast path leaves to the reference's code (an empty GT
    subfield, a multi-character allele, a space inside a later field)."""
    rng = random.Random(9)
    formats = {
        'GT:DS': lambda g: g + ':' + rng.choice(['0.5', '1.998', '0', '.']),
        'GT:AD:DP:GQ': lambda g: '%s:%d,%d:%d:%d' % (g, rng.randint(0, 30), rng.randint(0, 30), rng.randint(0, 60),
                                                     rng.randint(0, 99)),
        'GT:': lambda g: g + ':',
        'GT:DS:GP': lambda g: rng.choice([g, g + ':0.1', g + ':0.1:0.9,0.1,0']),        # trailing fields dropped
    }
    records = []
    for i in range(60):
        fmt = rng.choice(list(formats))
        alt = rng.choice(['G', 'GT', 'C,G', 'A,C,T'])
        n_alt = alt.count(',') + 1
        tokens = [formats[fmt](rng.choice(['0|0', '0|0', '0|1', '1|0', '%d|%d' % (rng.randint(0, n_alt), n_alt),
                                           '1', '0', '.', '.|.', '0/1', './.']))
                  for _ in SAMPLES]
        if i % 7 == 0:
            tokens[rng.randrange(len(tokens))] = '.'                                    # the whole token missing
        if i % 11 == 0:
            tokens[rng.randrange(len(tokens))] = rng.choice([':0.5', '10|0:0.5', '0|1:a b'])
        afs = ','.join(str(round(rng.random() / 5, 3)) for _ in range(n_alt))
        records.append((100 + i * 10, rng.choice(['.', 'rs%d' % i]), 'A', alt, 'AC=1;AF=' + afs, fmt, tokens))
    records.append((900, '.', 'A', 'G', 'AF=0.1', 'GT:DS', ['0|0:0'] * 10))                # all REF
    records.append((910, '.', 'A', 'G', 'AF=0.1', 'GT:DS', ['0|1'] * 10))                  # no DS anywhere
    write(tmp_path / 'in.vcf.gz', vcf(SAMPLES, records))
    out = check(run_both(tmp_path, tmp_path / 'in.vcf.gz', af))
    assert 'GT:' not in out and '\tGT\t' in out


@pytest.mark.parametrize('case', ['extra_format', 'format_gt_second', 'spaces', 'multichar_allele', 'ten_alleles',
                                  'trailing_tab', 'double_tab', 'odd_separator', 'missing_genotypes'])
def test_lines_for_the_reference_code(tmp_path, case):
    """Lines outside the fast path go through the reference's own code for that line; no fallback. ('extra_format',
    GT:DS, is on the fast path since GT comes first; kept here as a multi-allelic line among GT ones.)"""
    rng = random.Random(4)
    good = [(100 + i, '.', 'A', 'G', 'AF=0.2', 'GT', random_genotypes(rng, 10)) for i in range(3)]
    text = vcf(SAMPLES, good)
    gts = random_genotypes(rng, 10)
    line = '\t'.join(['22', '300', '.', 'A', 'G,T', '.', 'PASS', 'AF=0.1,0.3', 'GT'] + gts)
    if case == 'extra_format':
        line = '\t'.join(['22', '300', '.', 'A', 'G,T', '.', 'PASS', 'AF=0.1,0.3', 'GT:DS']
                         + [g + ':0.5' for g in gts])
    elif case == 'format_gt_second':
        line = '\t'.join(['22', '300', '.', 'A', 'G', '.', 'PASS', 'AF=0.1', 'DS:GT'] + ['0.5:' + g for g in gts])
    elif case == 'spaces':
        line = line.replace('\t', ' ', 12) + '  '          # the reference splits on any whitespace
    elif case == 'multichar_allele':
        line = line.replace(gts[3], '10|0', 1) if gts[3] in line else line + '\t10|0'
    elif case == 'ten_alleles':
        line = '\t'.join(['22', '300', '.', 'A', ','.join(['C'] * 10), '.', 'PASS', 'AF=' + ','.join(['0.1'] * 10),
                          'GT'] + ['%d|0' % rng.randint(0, 9) for _ in range(10)])
    elif case == 'trailing_tab':
        line = line + '\t'
    elif case == 'double_tab':
        line = line.replace('\tGT\t', '\tGT\t\t')
    elif case == 'odd_separator':
        line = line[:-3] + '0-1'
    elif case == 'missing_genotypes':
        line = '\t'.join(['22', '300', '.', 'A', 'G', '.', 'PASS', 'AF=0.1', 'GT'])
    write(tmp_path / 'in.vcf', text + line + '\n' + vcf(SAMPLES, good).split('\n', 3)[3])
    check(run_both(tmp_path, tmp_path / 'in.vcf', '0.05'))


@pytest.mark.parametrize('case', ['carriage_return', 'non_ascii_header', 'non_ascii_line'])
def test_whole_file_to_the_reference(tmp_path, case):
    rng = random.Random(5)
    records = [(100 + i, '.', 'A', 'G', 'AF=0.2', 'GT', random_genotypes(rng, 10)) for i in range(5)]
    text = vcf(SAMPLES, records)
    if case == 'carriage_return':
        text = text.replace('\n', '\r\n')
    elif case == 'non_ascii_header':
        text = text.replace('##contig=<ID=22>', '##source=Mänttä')
    else:
        text = text.replace('AF=0.2', 'AF=0.2;NOTE=é', 1)
    write(tmp_path / 'in.vcf', text)
    check(run_both(tmp_path, tmp_path / 'in.vcf'), fallback=True)


def test_custom_af_field(tmp_path):
    rng = random.Random(6)
    records = [(100 + i, '.', 'A', rng.choice(['G', 'C,G']), 'AC=1;MYAF=%s;AF=0.9' % rng.choice(['0.01', '0.2', '0.01,0.3']),
                'GT', random_genotypes(rng, 10)) for i in range(20)]
    records = [r for r in records if (',' in r[3]) == (',' in r[4])]
    write(tmp_path / 'in.vcf', vcf(SAMPLES, records))
    check(run_both(tmp_path, tmp_path / 'in.vcf', '0.1', ['-af_field', 'MYAF']))


def test_both_stop_on_what_stops_the_reference(tmp_path):
    """A space inside INFO shifts the reference's fields, and it stops on FORMAT; so does src/."""
    rng = random.Random(7)
    records = [(100, '.', 'A', 'G', 'AF=0.2;NOTE=a b', 'GT', random_genotypes(rng, 10))]
    write(tmp_path / 'in.vcf', vcf(SAMPLES, records))
    results = run_both(tmp_path, tmp_path / 'in.vcf')
    assert results['legacy'][0] != 0 and results['fast'][0] != 0
    assert "'GT' is not in list" in results['legacy'][3] and "'GT' is not in list" in results['fast'][3]


# ------------------------------------ -o - (F3: prohap.py -filter_glob reads it through a pipe) ------------------------

def run_stdout(pattern, af='0', env=None):
    return subprocess.run([sys.executable, str(TREES['fast'] / 'vcf_filter_fix.py'), '-i', str(pattern), '-chr', '22',
                           '-af', af, '-o', '-'], capture_output=True,
                          env=dict(os.environ, PYTHONHASHSEED='0', **(env or {})))


def test_stdout_is_the_file_uncompressed(tmp_path):
    rng = random.Random(8)
    def records(start):
        out = []
        for i in range(30):
            alt, af = rng.choice([('T', 'AF=0.2'), ('T', 'AF=0.01'), ('G,T', 'AF=0.2,0.01')])
            out.append((start + i, '.', 'C', alt, af, 'GT',
                        random_genotypes(rng, 10, ('0|0', '0|1', '2|1', '1', '.|0'))))
        return out
    write(tmp_path / 'c.1.vcf.gz', vcf(SAMPLES, records(100)))
    write(tmp_path / 'c.2.vcf', vcf(SAMPLES, records(300), header=HEADER + ['##second=1']).rstrip('\n'))
    write(tmp_path / 'c.3.vcf.gz', '')
    results = run_both(tmp_path, tmp_path / 'c.*', '0.05')
    expected = check(results)
    proc = run_stdout(tmp_path / 'c.*', '0.05')
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.decode() == expected
    # the messages go to stderr, nothing else does
    assert proc.stderr.decode().splitlines() == results['fast'][1].splitlines()
    assert not (tmp_path / '-').exists()


@pytest.mark.parametrize('fallback', [None, 'error'])
def test_stdout_stops_instead_of_falling_back(tmp_path, fallback):
    rng = random.Random(9)
    records = [(100 + i, '.', 'A', 'G', 'AF=0.2', 'GT', random_genotypes(rng, 10)) for i in range(5)]
    write(tmp_path / 'in.vcf', vcf(SAMPLES, records).replace('\n', '\r\n'))
    proc = run_stdout(tmp_path / 'in.vcf', env={'PROHAP_FALLBACK': fallback} if fallback else None)
    assert proc.returncode == 3
    assert b'the filter stops on input it does not handle' in proc.stderr and FALLBACK.encode() not in proc.stderr
    assert not (tmp_path / '-').exists() and not (Path.cwd() / '-').exists()
