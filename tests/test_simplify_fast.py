"""src/fasta_simplify_headers.py (S1) against legacy/src/fasta_simplify_headers.py, as scripts, on crafted inputs:

    pytest tests/test_simplify_fast.py

The reference reads the whole final FASTA into a dict, the haplotype table into pandas, and writes the header table
with one to_csv. The fast script streams the FASTA, keeps only HaplotypeID -> TranscriptID, and writes the header
table one block of records at a time. Each case runs both scripts on the same FASTA, haplotype table and a small
gffutils DB, and checks that the fast one ran without handing over to the reference when it should, handed over when
it should, and that both end the same way with the same bytes in both outputs.
"""
import bz2
import concurrent.futures
import gzip
import os
import subprocess
import sys
from pathlib import Path

import pytest

gffutils = pytest.importorskip('gffutils')

ROOT = Path(__file__).resolve().parent.parent
TREES = {'legacy': ROOT / 'legacy' / 'src', 'fast': ROOT / 'src'}
FALLBACK = 'ProHap-fast: falling back to the reference code'
BLOCK = 1 << 14             # src/fasta_simplify_headers.py: records per block
TABLE_CHUNK = 1 << 20       # rows per array of IDs while reading the table
COLUMNS = ['TranscriptID', 'chromosome', 'transcript_biotype', 'HaplotypeID', 'VCF_IDs', 'DNA_changes',
           'allele_frequencies', 'cDNA_changes', 'all_protein_changes', 'variant_types', 'protein_changes',
           'reading_frame', 'protein_prefix_length', 'start_missing', 'start_lost', 'splice_sites_affected',
           'occurrence_count', 'frequency', 'frequency_population', 'frequency_superpopulation']

# transcript: (gene_id or None, gene_name or None)
TRANSCRIPTS = {
    'ENST00000000001': ('ENSG00000000001', 'ALPHA'),
    'ENST00000000002': ('ENSG00000000002', None),           # no gene name: the gene ID is used
    'ENST00000000003': ('ENSG00000000001', 'ALPHA'),        # the same gene as ...01
    'ENST00000000004': ('ENSG00000000004', 'BETA'),
    'ENST00000000005': (None, None),                        # no gene ID: no gene at all
    'ENST00000000006': ('ENSG00000000006', 'GAMMA'),
}
TR = sorted(TRANSCRIPTS)


def gtf():
    lines = []
    for n, (tr, (gene, name)) in enumerate(sorted(TRANSCRIPTS.items())):
        attrs = ''
        if gene:
            attrs += 'gene_id "%s"; ' % gene
        attrs += 'transcript_id "%s";' % tr
        if name:
            attrs += ' gene_name "%s";' % name
        start = 1000 * (n + 1)
        lines.append('\t'.join(['21', 'ensembl', 'transcript', str(start), str(start + 500), '.', '+', '.', attrs]))
    return '\n'.join(lines) + '\n'


def record(tag, accession, ids, seq='MKVLA', rf='0'):
    return '>%s|%s|position_within_protein:0 start:0 matching_proteins:%s reading_frame:%s\n%s\n' % (
        tag, accession, ids, rf, seq)


def table_text(rows, columns=COLUMNS):
    """rows: (TranscriptID, HaplotypeID) or full rows."""
    out = ['\t'.join(columns)]
    for row in rows:
        if len(row) == 2:
            tr, hap = row
            full = [tr, '21', 'protein_coding', hap, 'rs1', '100:A>G', '0.01', '10:A>G', '4:V>4:D', 'SAV',
                    '4:V>4:D', '0', '0', 'False', 'False', '-', '3', '3e-06', 'GBR:0.001', 'EUR:0.001']
            row = [full[COLUMNS.index(c)] if c in COLUMNS else 'x' for c in columns]
        out.append('\t'.join(row))
    return '\n'.join(out) + '\n'


def haplotypes(n, start=0, prefix='haplo_chr21_'):
    """n (transcript, haplotype ID) rows."""
    return [(TR[(start + i) % 4 if (start + i) % 7 else 3], '%s%x' % (prefix, start + i)) for i in range(n)]


def basic_fasta(n_hap=40, prefix='haplo_chr21_', var_prefix='var_'):
    parts = [record('generic_ensref', 'prot_0', 'ENST00000000001;ENST00000000003'),
             record('generic_ensref', 'prot_1', 'ENST00000000002,ENST00000000004'),
             record('generic_ensref', 'prot_2', 'ENST00000000005'),
             record('generic_cont', 'CONT_1', 'P12345;ENST00000000001,haplo_x;Q99999'),
             record('generic_var', 'prot_3', '%sENST00000000006_1:A>G;%sENST00000000004_7:C>T' % (var_prefix,
                                                                                                   var_prefix)),
             record('generic_other', 'prot_4', 'whatever'),
             '>generic_ensref|prot_5|position_within_protein:0 start:0 matching_proteins:ENST00000000006 '
             'reading_frame:-\nMKV\nLAA\n\nWWW\n',
             'ensref|prot_6 position_within_protein:3 start:1 matching_proteins:ENST00000000001 reading_frame:1\n'
             'MMM\n']
    parts[-1] = '>' + parts[-1]
    for i in range(n_hap):
        ids = ['%s%x' % (prefix, i)] + (['%s%x' % (prefix, i + 1)] if i % 3 == 0 else [])
        parts.append(record('generic_enshap', 'prot_h%d' % i, ';'.join(ids) if i % 2 else ','.join(ids),
                            rf=','.join('0' * len(ids))))
    return ''.join(parts)


def big_fasta(n):
    """n haplotype records, in blocks of BLOCK."""
    parts = [record('generic_ensref', 'prot_r', 'ENST00000000001')]
    for i in range(n):
        parts.append(record('generic_enshap', 'p%x' % i, 'haplo_chr21_%x' % i, seq='MK' + 'A' * (i % 50)))
    return ''.join(parts)


def chunk_edge_table():
    """More rows than TABLE_CHUNK, with the longest IDs after the first chunk; the FASTA looks up IDs on both sides."""
    rows = haplotypes(TABLE_CHUNK + 50)
    rows[-3] = ('ENST00000000006', 'haplo_chr21_a_much_longer_identifier_than_the_others')
    return table_text(rows, ['HaplotypeID', 'TranscriptID', 'x'])


def chunk_edge_fasta():
    ids = ['haplo_chr21_%x' % i for i in (0, TABLE_CHUNK - 1, TABLE_CHUNK, TABLE_CHUNK + 49)]
    ids.append('haplo_chr21_a_much_longer_identifier_than_the_others')
    return ''.join(record('generic_enshap', 'p%d' % i, x) for i, x in enumerate(ids))


def replace_line(text, n, line):
    lines = text.split('\n')
    lines[n] = line
    return '\n'.join(lines)


T = table_text(haplotypes(60))
F = basic_fasta()

# name: (fasta, table or None, extra arguments, fasta suffix, output suffix, header suffix, fast?, does it stop?)
CASES = {
    'gzipped': (F, T, [], '.fa.gz', '.fasta.gz', '.tsv.gz', True, False),
    'plain': (F, T, [], '.fa', '.fasta', '.tsv', True, False),
    'block_edges': (big_fasta(2 * BLOCK + 3), table_text(haplotypes(2 * BLOCK + 3)), [], '.fa.gz', '.fasta.gz',
                    '.tsv.gz', True, False),
    'exactly_one_block': (big_fasta(BLOCK - 1), table_text(haplotypes(BLOCK)), [], '.fa', '.fasta', '.tsv', True,
                          False),
    'table_chunk_edge': (chunk_edge_fasta(), chunk_edge_table(), [], '.fa', '.fasta', '.tsv.gz', True, False),
    'no_table': (''.join(basic_fasta(0)), None, [], '.fa', '.fasta', '.tsv', True, False),
    'no_haplotype_records': (basic_fasta(0), table_text([('ENST00000000001', 'NA'), ('12', '1e5'), ('x', 'NA')]),
                             [], '.fa', '.fasta', '.tsv', True, False),
    'empty_fasta': ('', T, [], '.fa', '.fasta', '.tsv', True, False),
    'header_only_table': (basic_fasta(0), table_text([]), [], '.fa', '.fasta', '.tsv', True, False),
    'ambiguous_other_columns': (F, table_text([(tr, hap, ) for tr, hap in haplotypes(60)]).replace(
        '\trs1\t', '\t2e5\t', 7).replace('\t3e-06\t', '\tTRUE\t', 5), [], '.fa', '.fasta', '.tsv', True, False),
    'unused_rows_numeric_transcript': (F, table_text(haplotypes(60) + [('12', 'haplo_unused')]), [], '.fa',
                                       '.fasta', '.tsv', True, False),
    'quotes_tabs_spaces': (F + '>generic_enshap|q"u\tote|position_within_protein:0 start:0 '
                           'matching_proteins:haplo_chr21_1 reading_frame:0\nMK\n'
                           '>generic_ens"ref| sp|position_within_protein:0 start:0 matching_proteins:ENST00000000004 '
                           'reading_frame:0\nMK\n', T, [], '.fa', '.fasta', '.tsv', True, False),
    'last_line_without_newline': (F + '>generic_ensref|prot_z|position_within_protein:0 start:0 '
                                  'matching_proteins:ENST00000000001 reading_frame:-\nMKVL', T, [], '.fa', '.fasta',
                                  '.tsv', True, False),
    'custom_prefixes': (basic_fasta(prefix='hp_', var_prefix='v_'), table_text(haplotypes(60, prefix='hp_')),
                        ['-hap_prefix', 'hp_', '-var_prefix', 'v_'], '.fa', '.fasta', '.tsv', True, False),
    'output_not_gz_suffix': (F, T, [], '.fa', '.fasta.bz2', '.tsv', True, False),
    # handed to the reference, which writes the outputs
    'duplicate_accession': (F + record('generic_ensref', 'prot_1', 'ENST00000000006'), T, [], '.fa', '.fasta', '.tsv',
                            False, False),
    'carriage_return': (F.replace('\n', '\r\n'), T, [], '.fa', '.fasta', '.tsv', False, False),
    'non_ascii_fasta': (F.replace('prot_2', 'prot_é'), T, [], '.fa', '.fasta', '.tsv', False, False),
    'quote_in_table': (F, T.replace('\trs1\t', '\t"rs1"\t', 1), [], '.fa', '.fasta', '.tsv', False, False),
    'short_table_line': (F, replace_line(T, 55, 'ENST00000000001\t21'), [], '.fa', '.fasta', '.tsv', False, False),
    'blank_table_line': (F, replace_line(T, 50, ''), [], '.fa', '.fasta', '.tsv', False, False),
    'na_haplotype_id': (F, table_text(haplotypes(60) + [('ENST00000000001', 'NA')]), [], '.fa', '.fasta', '.tsv',
                        False, False),
    'numeric_haplotype_id': (F, table_text(haplotypes(60) + [('ENST00000000001', '1e5')]), [], '.fa', '.fasta',
                             '.tsv', False, False),
    'bz2_header': (F, T, [], '.fa', '.fasta', '.tsv.bz2', False, False),
    # the reference stops; so does the fast script, through the reference
    'missing_haplotype': (F, table_text(haplotypes(30)), [], '.fa', '.fasta', '.tsv', False, True),
    'haplotypes_without_table': (F, None, [], '.fa', '.fasta', '.tsv', False, True),
    'duplicate_haplotype_id': (F, table_text(haplotypes(60) + [('ENST00000000004', 'haplo_chr21_5')]), [], '.fa',
                               '.fasta', '.tsv', False, True),
    'numeric_transcript_id': (F, table_text([('12' if hap == 'haplo_chr21_4' else tr, hap)
                                             for tr, hap in haplotypes(60)]), [], '.fa', '.fasta', '.tsv', False,
                              True),
    'transcript_not_in_db': (F + record('generic_ensref', 'prot_x', 'ENST00000000099'), T, [], '.fa', '.fasta',
                             '.tsv', False, True),
    'no_matching_proteins': (F + '>generic_ensref|prot_y|position_within_protein:0 start:0\nMK\n', T, [], '.fa',
                             '.fasta', '.tsv', False, True),
    'variant_without_enst': (F + record('generic_var', 'prot_v', 'var_X_1'), T, [], '.fa', '.fasta', '.tsv', False,
                             True),
    'no_haplotype_id_column': (F, T.replace('HaplotypeID', 'HaplotypeName', 1), [], '.fa', '.fasta', '.tsv', False,
                               True),
    'no_transcript_id_column': (F, T.replace('TranscriptID', 'Transcript', 1), [], '.fa', '.fasta', '.tsv', False,
                                True),
}


def write(path, text):
    data = text.encode('utf-8')
    if str(path).endswith('.gz'):
        data = gzip.compress(data)
    elif str(path).endswith('.bz2'):
        data = bz2.compress(data)
    path.write_bytes(data)


def read(path):
    data = path.read_bytes()
    if data[:2] == b'\x1f\x8b':
        return gzip.decompress(data)
    if data[:3] == b'BZh':
        return bz2.decompress(data)
    return data


def run_case(case, db, name):
    fasta, table, extra, suffix, out_suffix, header_suffix, _, _ = CASES[name]
    fasta_path = case / ('final' + suffix)
    write(fasta_path, fasta)
    command = ['-i', str(fasta_path), '-db', str(db)] + extra
    if table is not None:
        table_path = case / 'haplo.tsv.gz'
        write(table_path, table)
        command += ['-hap_tsv', str(table_path)]
    env = dict(os.environ, PYTHONHASHSEED='0', PYTHONNOUSERSITE='1')
    env.pop('PROHAP_FALLBACK', None)
    result = {}
    for arm, src in TREES.items():
        out = case / arm
        out.mkdir()
        outputs = [out / ('simplified' + out_suffix), out / ('header' + header_suffix)]
        proc = subprocess.run([sys.executable, str(src / 'fasta_simplify_headers.py')] + command +
                              ['-o', str(outputs[0]), '-header', str(outputs[1])],
                              capture_output=True, text=True, env=env)
        files = [read(p) if p.exists() else None for p in outputs] if proc.returncode == 0 else None
        result[arm] = proc, files
    return result


@pytest.fixture(scope='module')
def results(tmp_path_factory):
    """Every case, run by both trees, eight cases at a time."""
    root = tmp_path_factory.mktemp('simplify')
    db = root / 'annotation.db'
    gffutils.create_db(gtf(), dbfn=str(db), from_string=True, force=True, keep_order=True, merge_strategy='merge',
                       sort_attribute_values=True, disable_infer_genes=True, disable_infer_transcripts=True)

    def one(name):
        case = root / name
        case.mkdir()
        return name, run_case(case, db, name)
    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        return dict(pool.map(one, CASES))


@pytest.mark.parametrize('name', list(CASES))
def test_case(results, name):
    fast_expected, stops = CASES[name][6:8]
    legacy, legacy_files = results[name]['legacy']
    fast, files = results[name]['fast']
    assert (FALLBACK in fast.stderr) == (not fast_expected), fast.stderr
    if stops:
        assert legacy.returncode != 0 and fast.returncode != 0, (legacy.stderr, fast.stderr)
    else:
        assert legacy.returncode == 0, legacy.stderr
        assert fast.returncode == 0, fast.stderr
        assert None not in files
        assert files == legacy_files


def test_cases_cover_the_outputs(results):
    """The fast cases write records of every kind the reference distinguishes."""
    _, files = results['gzipped']['fast']
    fasta = files[0].decode()
    for text in (' GN=ALPHA\n', ' GN=ENSG00000000002\n', ' GN=missing\n', ' CONTAMINANT GN=P12345;Q99999\n',
                 '.0 GN=', '.1 GN=', ' GN=GAMMA\n', 'MKVLAAWWW\n'):
        assert text in fasta
    _, files = results['quotes_tabs_spaces']['fast']
    assert b'"q""u\tote"' in files[1]
