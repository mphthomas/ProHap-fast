"""prohap.py -fix_bugs 1 (config fix_bugs: true): the fixes of the reference's haplotype bugs (IMPLEMENTATION_NOTES
section 9), each against hand-worked expectations. With -fix_bugs 0 the fast code is checked against the reference
byte for byte elsewhere (test_prohap_fast.py, tests/legacy_vs_fast.sh, verify); here too, on the crafted chromosome.

The crafted chromosome has two transcripts with a designed cDNA each, both with a 6-base 5' UTR, 20 codons and a stop,
split over two exons:

    TPLUS   + strand, exons 101-140 and 181-233; cDNA index i is at 101 + i (i < 40) or 181 + i - 40
            GCCACC ATG GCT AAA CTG GTT CAG GAA TGC CGT TTC GAT AGC CCA TGG GAG CAC ATC TAC AAC GGA TAA UTR
                   M   A   K   L   V   Q   E   C   R   F   D   S   P   W   E   H   I   Y   N   G
    TMINUS  - strand, exons 301-353 and 394-433; cDNA index i is at 433 - i (i < 40) or 353 - (i - 40), complemented
            GCCACC ATG TCC GAA CGC AAG TTA CCT GGC GAC TAT CAA ACG GTG TGG AGA CTC GCA ATT TTC CAG TGA UTR
                   M   S   E   R   K   L   P   G   D   Y   Q   T   V   W   R   L   A   I   F   Q
    codon j covers cDNA 6 + 3j .. 8 + 3j; the 3' UTR is GGGCCCAAATTTGGGCCCAAATTT, without an ATG

Each case is the first copy of one sample (the second copy is the reference), carrying the case's VCF records (genome
coordinates, REF from the genome). The expected outcome of each case is written by hand for each set of fixes: a
protein (ignore_UTR 1: from the start codon to the first stop), 'removed' (the copy is in the removed-samples log,
stage 1) or 'ref' (no row: the copy is coded as the reference).

    pytest tests/test_fix_bugs.py
"""
import gzip
import os
import random
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT / 'src') not in sys.path:
    sys.path.insert(0, str(ROOT / 'src'))
from modules import fixes  # noqa: E402
from modules.haplotypes import conflicting, trimmed_span  # noqa: E402

gffutils = pytest.importorskip('gffutils')

UTR3 = 'GGGCCCAAATTTGGGCCCAAATTT'
CDNA = {
    'TPLUS': 'GCCACC' + 'ATGGCTAAACTGGTTCAGGAATGCCGTTTCGATAGCCCATGGGAGCACATCTACAACGGATAA' + UTR3,
    'TMINUS': 'GCCACC' + 'ATGTCCGAACGCAAGTTACCTGGCGACTATCAAACGGTGTGGAGACTCGCAATTTTCCAGTGA' + UTR3,
}
REFERENCE_PROTEIN = {'TPLUS': 'MAKLVQECRFDSPWEHIYNG', 'TMINUS': 'MSERKLPGDYQTVWRLAIFQ'}
COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def genome():
    rng = random.Random(7)
    g = [rng.choice('ACGT') for _ in range(500)]
    plus, minus = CDNA['TPLUS'], CDNA['TMINUS']
    for i in range(93):
        g[(101 + i if i < 40 else 181 + i - 40) - 1] = plus[i]
        g[(433 - i if i < 40 else 353 - (i - 40)) - 1] = minus[i].translate(COMPLEMENT)
    return ''.join(g)


GENOME = genome()


def G(pos, n=1):
    return GENOME[pos - 1:pos - 1 + n]


def snv(pos, alt, ref):
    assert G(pos) == ref, (pos, G(pos), ref)        # the hand-worked base
    return (pos, ref, alt)


def deletion(anchor, n, ref):
    assert G(anchor, n + 1) == ref, (anchor, G(anchor, n + 1), ref)
    return (anchor, ref, ref[0])


def insertion(anchor, inserted, ref):
    assert G(anchor) == ref, (anchor, G(anchor), ref)
    return (anchor, ref, ref + inserted)


# case: (transcript, records, {fix name or 'off': expected from that fix on}). The fixes come in the order of
# fixes.ALL; a case keeps its expectation until a later fix names another.
CASES = {
    # 1: adjacent SNVs. + strand: two SNVs in one codon (GCT -> CAT, A2H); - strand: across codons 3 and 4 (E3D, R4C)
    'P1': ('TPLUS', [snv(110, 'C', 'G'), snv(111, 'A', 'C')],
           {'off': 'removed', 'conflicts': 'MHKLVQECRFDSPWEHIYNG'}),
    'M1': ('TMINUS', [snv(418, 'A', 'G'), snv(419, 'A', 'T')],
           {'off': 'removed', 'conflicts': 'MSDCKLPGDYQTVWRLAIFQ'}),
    # 2: a variant right after a deletion. + : L4 deleted, V5I; - : L6 deleted, K5N
    'P2': ('TPLUS', [deletion(115, 3, 'ACTG'), snv(119, 'A', 'G')],
           {'off': 'removed', 'conflicts': 'MAKLIQECRFDSPWEHIYNG',            # the + strand deletion: case 9
            'liftover': 'MAKIQECRFDSPWEHIYNG'}),
    'M2': ('TMINUS', [deletion(409, 3, 'GTAA'), snv(413, 'G', 'C')],
           {'off': 'removed', 'conflicts': 'MSERNPGDYQTVWRLAIFQ'}),
    # 3: a variant right before an indel's anchor base. + : C8Y, R9 deleted; - : D9 deleted, Y10C
    'P3': ('TPLUS', [snv(129, 'A', 'G'), deletion(130, 3, 'CCGT')],
           {'off': 'removed', 'conflicts': 'MAKLVQEYRFDSPWEHIYNG', 'liftover': 'MAKLVQEYFDSPWEHIYNG'}),
    'M3': ('TMINUS', [snv(399, 'C', 'T'), deletion(400, 3, 'AGTC')],
           {'off': 'removed', 'conflicts': 'MSERKLPGCQTVWRLAIFQ'}),
    # 4: an insertion followed by a SNV at the next base. + : K inserted after S12, P13S; - : W14C, P inserted after
    # it (in the genome the insertion comes first, on the cDNA the SNV)
    'P4': ('TPLUS', [insertion(182, 'AAA', 'C'), snv(183, 'T', 'C')],
           {'off': 'removed', 'conflicts': 'MAKLVQECRFDSKSWEHIYNG'}),
    'M4': ('TMINUS', [insertion(345, 'GGG', 'T'), snv(346, 'G', 'C')],
           {'off': 'removed', 'conflicts': 'MSERKLPGDYQTVCPRLAIFQ'}),
    # 5: a SNV on an indel's anchor base, which the trimmed forms combine. + : H16Q and I17 deleted; F10L and D
    # inserted after it; - : I18 deleted and F19V
    'P5': ('TPLUS', [snv(194, 'G', 'C'), deletion(194, 3, 'CATC')],
           {'off': 'removed', 'conflicts': 'MAKLVQECRFDSPWEQIYNG', 'liftover': 'MAKLVQECRFDSPWEQYNG'}),
    'P5b': ('TPLUS', [snv(136, 'A', 'C'), insertion(136, 'GAT', 'C')],
            {'off': 'removed', 'conflicts': 'MAKLVQECRLDSPWEHIYNG', 'liftover': 'MAKLVQECRLDDSPWEHIYNG'}),
    'M5': ('TMINUS', [snv(333, 'C', 'A'), deletion(333, 3, 'AAAT')],
           {'off': 'removed', 'conflicts': 'MSERKLPGDYQTVWRLAVQ'}),
    # 6-8 and two identical records stay conflicts
    'P6': ('TPLUS', [snv(120, 'C', 'T'), snv(120, 'A', 'T')], {'off': 'removed'}),
    'M6': ('TMINUS', [snv(405, 'A', 'C'), snv(405, 'T', 'C')], {'off': 'removed'}),
    'P7': ('TPLUS', [deletion(185, 3, 'ATGG'), snv(186, 'C', 'T')], {'off': 'removed'}),
    'M7': ('TMINUS', [deletion(400, 3, 'AGTC'), snv(401, 'A', 'G')], {'off': 'removed'}),
    'P8': ('TPLUS', [deletion(185, 3, 'ATGG'), deletion(186, 3, 'TGGG')], {'off': 'removed'}),
    'M8': ('TMINUS', [deletion(409, 3, 'GTAA'), deletion(410, 3, 'TAAC')], {'off': 'removed'}),
    'Pdup': ('TPLUS', [snv(125, 'C', 'G'), snv(125, 'C', 'G')], {'off': 'removed'}),
    # 9: indels alone. + : E15 deleted; G inserted after L4 (frameshift; the inserted G repeats the next base);
    # - : one T of TT deleted in F19 (frameshift; the deleted base repeats the anchor in the genome)
    'P9a': ('TPLUS', [deletion(188, 3, 'GGAG')], {'off': 'ref', 'liftover': 'MAKLVQECRFDSPWHIYNG'}),
    'P9b': ('TPLUS', [insertion(118, 'G', 'G')], {'off': 'ref', 'liftover': 'MAKLGSGMPFR'}),
    'M9': ('TMINUS', [deletion(332, 1, 'AA')], {'off': 'ref', 'liftover': 'MSERKLPGDYQTVWRLAISSEGPNLGPN'}),
    # 9: an allele that became the reference in a liftover (the genome has G, the record says A>G): still left out
    'P9c': ('TPLUS', [(110, 'A', 'G')], {'off': 'ref'}),
    # 10: stop gained: E7*; with a 5' UTR variant; with F10L after the new stop. - : Q11*
    'P10a': ('TPLUS', [snv(125, 'T', 'G')], {'off': 'ref', 'stop': 'MAKLVQ'}),
    'P10b': ('TPLUS', [snv(103, 'T', 'C'), snv(125, 'T', 'G')], {'off': 'ref', 'stop': 'MAKLVQ'}),
    'P10c': ('TPLUS', [snv(125, 'T', 'G'), snv(134, 'C', 'T')], {'off': 'ref', 'stop': 'MAKLVQ'}),
    'M10': ('TMINUS', [snv(397, 'A', 'G')], {'off': 'ref', 'stop': 'MSERKLPGDY'}),
    # 10 by a frameshift: TA inserted after N19 makes the next codon TAG; the change reads '18:N>18:N*(+fs)', whose
    # ALT the reference measured as 7 residues long
    'P10d': ('TPLUS', [insertion(203, 'TA', 'C')], {'off': 'ref', 'stop': 'MAKLVQECRFDSPWEHIYN'}),
    # ... but not a frameshift right after the stop codon (TT inserted after TAA; the change reads '20:*>20:*L(+fs)')
    'P10e': ('TPLUS', [insertion(209, 'TT', 'A')], {'off': 'ref'}),
    # 11: start lost. M1L alone; M1L with L4M (an in-frame ATG); GGG inserted in the 5' UTR with M1I; AAA inserted
    # right after the ATG; - : M1T (no ATG downstream); one 5' UTR base deleted right before the ATG, with S2A
    # (ATG at cDNA 26 in frame 2: M P F R *; the new ATG at 15 in frame 0)
    'P11a': ('TPLUS', [snv(107, 'C', 'A')], {'off': 'ref', 'start': 'MPFR'}),
    'P11b': ('TPLUS', [snv(107, 'C', 'A'), snv(116, 'A', 'C')], {'off': 'ref', 'start': 'MVQECRFDSPWEHIYNG'}),
    'P11c': ('TPLUS', [insertion(103, 'GGG', 'C'), snv(109, 'A', 'G')],
             {'off': 'IAKLVQECRFDSPWEHIYNG', 'start': 'MPFR'}),        # the reference misses the lost start
    'P11d': ('TPLUS', [insertion(109, 'AAA', 'G')], {'off': 'ref', 'start': 'MKAKLVQECRFDSPWEHIYNG'}),
    'M11': ('TMINUS', [snv(426, 'G', 'A')], {'off': 'ref'}),
    'M11b': ('TMINUS', [deletion(427, 1, 'TG'), snv(424, 'C', 'A')], {'off': 'ref', 'start': 'MAERKLPGDYQTVWRLAIFQ'}),
    # a frameshift whose first codon translates the same: A inserted in AAA (K3); on the - strand T inserted before
    # TTC (F19 stays F as TTT). The reference calls them synonymous; until 'liftover' it also drops them
    'P12': ('TPLUS', [insertion(113, 'A', 'A')], {'off': 'ref', 'frameshift': 'MAKTGSGMPFR'}),
    'M12': ('TMINUS', [insertion(333, 'A', 'A')], {'off': 'ref', 'frameshift': 'MSERKLPGDYQTVWRLAIFPVRAQIWAQI'}),
}
NAMES = list(CASES)


def expected(case, applied):
    """The expectation of a case with the fixes `applied` (all of fixes.ALL up to some point)."""
    table = CASES[case][2]
    result = table['off']
    for name in fixes.ALL:
        if name in applied and name in table:
            result = table[name]
    return result


def write_inputs(path):
    """The crafted chromosome: gffutils DB, cDNA FASTA, transcript list, metadata and VCF."""
    gtf = []
    exons = {'TPLUS': ('+', [(101, 140), (181, 233)], (107, 109), (207, 209)),
             'TMINUS': ('-', [(301, 353), (394, 433)], (425, 427), (325, 327))}
    for tid, (strand, spans, start, stop) in exons.items():
        attrs = 'gene_id "G%s"; transcript_id "%s"; transcript_biotype "protein_coding";' % (tid, tid)
        row = lambda kind, s, e: '\t'.join(['22', 'test', kind, str(s), str(e), '.', strand, '.', attrs])  # noqa
        gtf.append(row('transcript', spans[0][0], spans[-1][1]))
        gtf.extend(row('exon', s, e) for s, e in spans)
        gtf.append(row('start_codon', *start))
        gtf.append(row('stop_codon', *stop))
    (path / 'in.gtf').write_text('\n'.join(gtf) + '\n')
    gffutils.create_db(str(path / 'in.gtf'), dbfn=str(path / 'in.db'), force=True, keep_order=True,
                       merge_strategy='merge', sort_attribute_values=True, disable_infer_genes=True,
                       disable_infer_transcripts=True)
    (path / 'cdna.fa').write_text(''.join('>%s.1 cdna chromosome:test\n%s\n' % (t, s) for t, s in CDNA.items()))
    (path / 'transcripts.csv').write_text('chromosome,transcriptID\n22,TPLUS\n22,TMINUS\n')
    samples = ['S%02d' % (i + 1) for i in range(len(NAMES))]
    (path / 'samples.tsv').write_text('Sample name\tSex\tPopulation code\tSuperpopulation code\n' +
                                      ''.join('%s\tfemale\tGBR\tEUR\n' % s for s in samples))
    lines = []
    for n, name in enumerate(NAMES):
        for k, (pos, ref, alt) in enumerate(CASES[name][1]):
            gts = ['1|0' if m == n else '0|0' for m in range(len(NAMES))]
            lines.append((pos, n, k, ['22', str(pos), '%s_%d' % (name, k), ref, alt, '.', 'PASS', 'AF=0.01', 'GT']
                          + gts))
    lines.sort(key=lambda x: x[:3])
    text = '##fileformat=VCFv4.2\n#' + '\t'.join(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
                                                  'FORMAT'] + samples) + '\n'
    text += ''.join('\t'.join(line[3]) + '\n' for line in lines)
    with gzip.open(path / 'in.vcf.gz', 'wt') as f:
        f.write(text)
    return {name: samples[n] + ':1' for n, name in enumerate(NAMES)}


def run(path, script, out, fix_bugs=None, fix_names=None, extra=()):
    out.mkdir()
    env = dict(os.environ, PYTHONHASHSEED='0', PYTHONNOUSERSITE='1', PROHAP_FALLBACK='error')
    env.pop('PROHAP_FIXES', None)
    if fix_names is not None:
        env['PROHAP_FIXES'] = ','.join(fix_names)
    command = [sys.executable, str(ROOT / 'tools' / 'run_fork.py'), str(script), '-i', str(path / 'in.vcf.gz'),
               '-db', str(path / 'in.db'), '-transcripts', str(path / 'transcripts.csv'), '-cdna',
               str(path / 'cdna.fa'), '-s', str(path / 'samples.tsv'), '-chr', '22', '-min_hap_freq', '-1',
               '-min_hap_count', '0', '-acc_prefix', 'enshap_22', '-id_prefix', 'haplo_chr22', '-threads', '1',
               '-log', 'prohap.log', '-tmp_dir', 'tmp', '-output_csv', 'haplo.tsv', '-output_fasta', 'haplo.fa',
               '-output_cdna_fasta', 'cdna.fa'] + list(extra)
    if fix_bugs is not None:
        command += ['-fix_bugs', str(fix_bugs)]
    (out / 'tmp').mkdir()
    proc = subprocess.run(command, cwd=out, capture_output=True, text=True, env=env)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return proc


def outcome(out):
    """{(transcript, copy): protein, 'removed'}, and the table rows by transcript."""
    sequences = {}
    lines = (out / 'haplo.fa').read_text().split('\n')
    for header, seq in zip(lines[0::2], lines[1::2]):
        for h in header.split('matching_proteins:')[1].split()[0].split(';'):
            sequences[h] = seq
    result = {}
    rows = []
    table = (out / 'haplo.tsv').read_text().split('\n')
    columns = table[0].split('\t')
    for line in table[1:]:
        if line:
            row = dict(zip(columns, line.split('\t')))
            rows.append(row)
            for copy in row['samples'].split(';'):
                result[(row['TranscriptID'], copy)] = sequences[row['HaplotypeID']]
    for line in (out / 'prohap.log').read_text().split('\n')[1:]:
        if line:
            transcript, _, removed = line.split('\t')
            for copy in removed.split(';') if removed else []:
                assert (transcript, copy) not in result
                result[(transcript, copy)] = 'removed'
    return result, rows


@pytest.fixture(scope='module')
def inputs(tmp_path_factory):
    path = tmp_path_factory.mktemp('crafted')
    return path, write_inputs(path)


def levels():
    """Every prefix of fixes.ALL: the fixes as they were added, step by step."""
    return [tuple(fixes.ALL[:k]) for k in range(len(fixes.ALL) + 1)]


@pytest.fixture(scope='module')
def results(inputs, tmp_path_factory):
    path, copies = inputs
    root = tmp_path_factory.mktemp('runs')
    runs = {}
    for level in levels():
        out = root / ('fixes-' + '-'.join(level or ('none',)))
        run(path, ROOT / 'src' / 'prohap.py', out, 1, level)
        runs[level] = outcome(out)
    runs['dir'] = out           # the run with every fix
    return copies, runs


@pytest.mark.parametrize('level', levels(), ids=lambda level: '+'.join(level) or 'none')
def test_cases(results, level):
    """Every case, with the fixes up to each step: the hand-worked outcome, and no other copy affected."""
    copies, runs = results
    got, _ = runs[level]
    wrong = {}
    for name in NAMES:
        transcript = CASES[name][0]
        want = expected(name, level)
        have = got.get((transcript, copies[name]), 'ref')
        if have != want:
            wrong[name] = (have, want)
    assert not wrong, wrong
    # only the cases' first copies, on their own transcript
    assert {k for k in got} <= {(CASES[n][0], copies[n]) for n in NAMES}


def test_no_fix_is_the_reference(inputs, tmp_path):
    """-fix_bugs 0 writes the reference's output byte for byte, and gives the 'off' expectations; so does
    -fix_bugs 1 with PROHAP_FIXES empty."""
    path, copies = inputs
    run(path, ROOT / 'legacy' / 'src' / 'prohap.py', tmp_path / 'legacy')
    run(path, ROOT / 'src' / 'prohap.py', tmp_path / 'off', 0)
    run(path, ROOT / 'src' / 'prohap.py', tmp_path / 'none', 1, ())
    for name in ('haplo.tsv', 'haplo.fa', 'cdna.fa', 'prohap.log'):
        legacy = (tmp_path / 'legacy' / name).read_bytes()
        assert (tmp_path / 'off' / name).read_bytes() == legacy, name
        assert (tmp_path / 'none' / name).read_bytes() == legacy, name
    got, _ = outcome(tmp_path / 'off')
    assert {n: got.get((CASES[n][0], copies[n]), 'ref') for n in NAMES} == {n: CASES[n][2]['off'] for n in NAMES}


def test_rows_keep_the_original_records(results):
    """Case 5: both records stay in Changes and VCF_IDs as the VCF gives them (the - strand in reverse order)."""
    copies, runs = results
    _, rows = runs[tuple(fixes.ALL)]
    m5 = next(r for r in rows if copies['M5'] in r['samples'].split(';'))
    assert m5['VCF_IDs'] == 'M5_1;M5_0' and m5['DNA_changes'] == '333:AAAT>A;333:A>C'


def test_stop_gain_rows(results):
    """Case 10: the stop-gain is in the row; a 5' UTR variant and a variant after the new stop are not."""
    copies, runs = results
    if 'stop' not in fixes.ALL:
        pytest.skip('no stop fix yet')
    _, rows = runs[tuple(fixes.ALL)]
    row = {name: next(r for r in rows if copies[name] in r['samples'].split(';')) for name in ('P10a', 'P10b', 'P10c')}
    assert [(r['DNA_changes'], r['VCF_IDs']) for r in row.values()] == \
        [('125:G>T', 'P10a_0'), ('125:G>T', 'P10b_1'), ('125:G>T', 'P10c_0')]
    assert all(r['protein_changes'] == '6:E>6:*' and r['variant_types'] == 'stop_gained' for r in row.values())


def test_rescanned_start_rows(results):
    """Case 11: a re-scanned start shows start_lost True and its reading frame in the table, start:0 and that frame
    in the FASTA (ignore_UTR 1: the sequence begins at the new start); the change that makes the new ATG stays in the
    row; a start that is not lost after all is not flagged."""
    copies, runs = results
    if 'start' not in fixes.ALL:
        pytest.skip('no start fix yet')
    _, rows = runs[tuple(fixes.ALL)]
    row = {name: next(r for r in rows if copies[name] in r['samples'].split(';'))
           for name in ('P11a', 'P11b', 'P11c', 'P11d', 'M11b')}
    assert [row[n]['start_lost'] for n in row] == ['True', 'True', 'True', 'False', 'False']
    assert [row[n]['reading_frame'] for n in row] == ['2', '0', '2', '0', '2']
    assert row['P11a']['variant_types'] == 'start_lost' and row['P11a']['protein_changes'].startswith('0:M>')
    assert row['P11b']['DNA_changes'] == '107:A>C;116:C>A'
    assert row['P11c']['DNA_changes'] == '109:G>A'      # the 5' UTR insertion is left out as UTR
    fasta = (Path(runs['dir']) / 'haplo.fa').read_text()
    header = next(line for line in fasta.split('\n') if line.startswith('>') and row['P11a']['HaplotypeID'] in
                  line.split('matching_proteins:')[1].split()[0].split(';'))
    assert header.split()[1:] == ['start:0', 'reading_frame:2;2']     # P11a and P11c: the same protein


# ------------------------------------ the conflict rule (fix 'conflicts') ------------------------------------

@pytest.mark.parametrize('ref, alt, span', [
    ('A', 'G', (100, 101)), ('ACG', 'A', (101, 103)), ('A', 'AT', (101, 101)), ('AC', 'GT', (100, 102)),
    ('ACGT', 'AT', (101, 103)), ('ACG', 'ATG', (101, 102)), ('AC', 'ATC', (101, 101)), ('A', 'A', None),
    ('-', 'A', (100, 100)), ('T', '-', (100, 101)), ('AAA', 'A', (101, 103)),
])
def test_trimmed_span(ref, alt, span):
    """Shared prefix first, then shared suffix: a deletion in a run is the bases after the anchor."""
    assert trimmed_span(100, ref, alt) == span


@pytest.mark.parametrize('records, conflict', [
    ([(100, 'A', 'G'), (101, 'C', 'T')], False),             # 1 adjacent SNVs
    ([(100, 'ACG', 'A'), (103, 'T', 'C')], False),           # 2 right after a deletion
    ([(99, 'T', 'G'), (100, 'ACG', 'A')], False),            # 3 right before the anchor
    ([(100, 'A', 'AT'), (101, 'C', 'G')], False),            # 4 insertion, then the next base
    ([(100, 'A', 'T'), (100, 'ACG', 'A')], False),           # 5 SNV on the anchor
    ([(100, 'A', 'T'), (100, 'A', 'AT')], False),            # 5 SNV on an insertion's anchor
    ([(100, 'ACG', 'A'), (100, 'A', 'AT')], False),          # an insertion at a deletion's start
    ([(100, 'ACG', 'A'), (102, 'G', 'GT')], False),          # ... at its end
    ([(100, 'A', 'G'), (100, 'A', 'T')], True),              # 6 two ALTs
    ([(100, 'ACG', 'A'), (101, 'C', 'T')], True),            # 7 SNV inside a deletion
    ([(100, 'ACG', 'A'), (101, 'CGT', 'C')], True),          # 8 overlapping deletions
    ([(100, 'A', 'G'), (100, 'A', 'G')], True),              # the same record twice
    ([(100, 'A', 'AT'), (100, 'A', 'AG')], True),            # two insertions at one point
    ([(100, 'ACGT', 'A'), (101, 'C', 'CA')], True),          # an insertion inside a deletion
    ([(100, 'A', 'A'), (100, 'A', 'G')], False),             # a record that changes nothing
])
def test_conflicting(records, conflict):
    assert conflicting([trimmed_span(*r) for r in records]) == conflict


def test_selected():
    assert fixes.selected(0, {}) == frozenset()
    assert fixes.selected(0, {'PROHAP_FIXES': ','.join(fixes.ALL)}) == frozenset()
    assert fixes.selected(1, {}) == frozenset(fixes.ALL)
    assert fixes.selected(1, {'PROHAP_FIXES': ''}) == frozenset()
    with pytest.raises(ValueError, match='unknown fixes'):
        fixes.selected(1, {'PROHAP_FIXES': 'no_such_fix'})
