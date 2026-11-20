"""verify's comparison step (verify/compare.py) and the list of files it compares (verify.expected_outputs).

    pytest tests/test_compare.py
"""
import gzip
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'verify'))
import compare  # noqa: E402
import verify  # noqa: E402

TEXT = b'header\nline 2\nline 3\n'


def plain(path, data):
    path.write_bytes(data)
    return str(path)


def gz(path, data, mtime=None):
    with open(path, 'wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=mtime) as fh:
        fh.write(data)
    return str(path)


def test_identical_plain_files(tmp_path):
    r = compare.compare('t.tsv', plain(tmp_path / 'a', TEXT), plain(tmp_path / 'b', TEXT))
    assert r.status == compare.IDENTICAL and r.ok
    assert r.sha_a == r.sha_b == hashlib.sha256(TEXT).hexdigest()
    assert r.difference() == ''


def test_gz_files_compare_decompressed(tmp_path):
    """Two gzip files of the same text differ in their header timestamp, and are IDENTICAL."""
    a, b = gz(tmp_path / 'a.gz', TEXT, mtime=1), gz(tmp_path / 'b.gz', TEXT, mtime=2)
    assert Path(a).read_bytes() != Path(b).read_bytes()
    r = compare.compare('t.tsv.gz', a, b)
    assert r.status == compare.IDENTICAL
    assert r.sha_a == hashlib.sha256(TEXT).hexdigest()
    assert compare.sha256(a) == r.sha_a


def test_gz_vs_plain_is_different(tmp_path):
    r = compare.compare('t.gz', gz(tmp_path / 'a.gz', TEXT), plain(tmp_path / 'b.gz', TEXT))
    assert r.status == compare.DIFFERENT and not r.ok
    assert r.sha_a == r.sha_b and r.line is None
    assert r.difference() == 'legacy is gzip, fast is not'
    r = compare.compare('t.gz', plain(tmp_path / 'c.gz', TEXT), gz(tmp_path / 'd.gz', TEXT))
    assert r.difference() == 'fast is gzip, legacy is not'


def test_gzip_is_recognised_by_content_not_name(tmp_path):
    assert compare.compare('t', gz(tmp_path / 'a.txt', TEXT), plain(tmp_path / 'b.txt', TEXT)).status == \
        compare.DIFFERENT
    r = compare.compare('t', plain(tmp_path / 'c.gz', TEXT), plain(tmp_path / 'd.gz', TEXT))
    assert r.status == compare.IDENTICAL


def test_first_differing_line(tmp_path):
    r = compare.compare('t', gz(tmp_path / 'a.gz', b'h\nx\t1\nsame\ny\n'),
                        gz(tmp_path / 'b.gz', b'h\nx\t1\nsame\nz\n'))
    assert r.status == compare.DIFFERENT
    assert (r.line, r.text_a, r.text_b) == (4, 'y', 'z')
    assert r.difference() == 'line 4'
    assert r.sha_a == hashlib.sha256(b'h\nx\t1\nsame\ny\n').hexdigest()
    assert r.sha_b == hashlib.sha256(b'h\nx\t1\nsame\nz\n').hexdigest()


def test_one_file_longer(tmp_path):
    r = compare.compare('t', plain(tmp_path / 'a', TEXT), plain(tmp_path / 'b', TEXT + b'extra\n'))
    assert (r.status, r.line, r.text_a, r.text_b) == (compare.DIFFERENT, 4, None, 'extra')
    r = compare.compare('t', plain(tmp_path / 'c', TEXT), plain(tmp_path / 'd', TEXT[:-1]))
    assert (r.status, r.line) == (compare.DIFFERENT, 3)


def test_empty_files(tmp_path):
    r = compare.compare('t', gz(tmp_path / 'a.gz', b''), gz(tmp_path / 'b.gz', b''))
    assert r.status == compare.IDENTICAL and r.sha_a == hashlib.sha256(b'').hexdigest()
    r = compare.compare('t', gz(tmp_path / 'c.gz', b''), gz(tmp_path / 'd.gz', b'x\n'))
    assert (r.status, r.line, r.text_a, r.text_b) == (compare.DIFFERENT, 1, None, 'x')


def test_long_lines_are_shortened(tmp_path):
    r = compare.compare('t', plain(tmp_path / 'a', b'A' * 500 + b'\n'), plain(tmp_path / 'b', b'A' * 499 + b'B\n'))
    assert r.line == 1 and r.text_a == 'A' * compare.SNIPPET + '...'


def test_missing_files(tmp_path):
    present = plain(tmp_path / 'a', TEXT)
    r = compare.compare('t', present, str(tmp_path / 'nope'))
    assert (r.status, r.note, r.sha_b) == (compare.MISSING, 'only in the legacy arm', None)
    assert r.sha_a == hashlib.sha256(TEXT).hexdigest() and not r.ok
    r = compare.compare('t', str(tmp_path / 'nope'), present)
    assert (r.status, r.note, r.sha_a) == (compare.MISSING, 'only in the fast arm', None)
    r = compare.compare('t', str(tmp_path / 'x'), str(tmp_path / 'y'))
    assert (r.status, r.note, r.ok) == (compare.MISSING, 'in neither arm', False)
    r = compare.compare('t', str(tmp_path / 'x'), str(tmp_path / 'y'), required=False)
    assert (r.status, r.ok) == (compare.ABSENT, True)


def test_corrupt_gzip(tmp_path):
    good = gz(tmp_path / 'a.gz', TEXT * 1000)
    data = Path(good).read_bytes()
    bad = plain(tmp_path / 'b.gz', data[:len(data) // 2])
    r = compare.compare('t', good, bad)
    assert r.status == compare.DIFFERENT and r.difference().startswith('cannot read')


# ------------------------------------ what is compared ------------------------------------

def base_config(**kw):
    cfg = {'ensembl_release': 113, 'chromosomes': ['21', 'X'], 'only_MANE_select': True,
           'working_dir_name_haplo': 'haplotypes_tmp', 'haplo_fasta_file': 'results/h.fa.gz',
           'haplo_table_file': 'results/h.tsv.gz', 'haplo_cdna_file': '', 'final_fasta_file': 'results/db.fa',
           'fasta_simplify_headers': True, 'use_ProVar': False, 'add_existing_haplo': False}
    cfg.update(kw)
    return cfg


def test_expected_outputs():
    out = verify.expected_outputs(base_config(), '/refs')
    assert out == [
        ('data/fasta/ensembl_reference_proteinDB_113_clean.fa.gz', True),
        ('data/fasta/crap_tagged.fa.gz', False),
        ('results/haplotypes_tmp/haplo_chr21.tsv.gz', True),
        ('results/haplotypes_tmp/haplo_chr21.fa.gz', True),
        ('log/prohap_chr21.log', False),
        ('results/haplotypes_tmp/haplo_chrX.tsv.gz', True),
        ('results/haplotypes_tmp/haplo_chrX.fa.gz', True),
        ('log/prohap_chrX.log', False),
        ('results/haplotypes_tmp/haplo_all.tsv.gz', True),
        ('results/h.tsv.gz', True),
        ('results/h_sampleIDs.tsv.gz', True),
        ('results/h.fa.gz', True),
        ('results/haplo_all_clean.fa.gz', True),
        ('results/ref_contam_vcf_haplo_all_clean.fa.gz', True),
        ('results/ref_contam_vcf_haplo_all_nodupl.fa.gz', True),
        ('results/db.fa', True),
        ('results/db_simplified.fasta', True),
        ('results/db_header.tsv', True),
    ]


def test_expected_outputs_options():
    names = lambda cfg: [n for n, _ in verify.expected_outputs(cfg, '/refs')]
    out = names(base_config(haplo_fasta_file='h.fa', haplo_table_file='out/h.v2.tsv', haplo_cdna_file='c.fa.gz',
                            final_fasta_file='out/db.v2.fa.gz', fasta_simplify_headers=False))
    assert 'results/haplotypes_tmp/haplo_chr21.fa' in out
    assert 'out/h.v2_sampleIDs.tsv' in out and 'c.fa.gz' in out
    assert out[-1] == 'out/db.v2.fa.gz'
    out = names(base_config(final_fasta_file='db.fa.gz'))
    assert out[-2:] == ['db_simplified.fasta.gz', 'db_header.tsv.gz']
    out = names(base_config(use_ProVar=True, var_table_file='v.tsv', var_fasta_file='v.fa', var_cdna_file='',
                            add_existing_haplo=True))
    assert ['v.tsv', 'v.fa', 'results/variants_all_clean.fa.gz', 'results/haplo_added_clean.fa.gz'] == \
        out[out.index('v.tsv'):out.index('v.tsv') + 4]
    # transcript lists built in the work directory are compared when they exist
    assert verify.expected_outputs(base_config(custom_transcript_list=''), '/refs')[0] == \
        ('data/included_transcripts.csv', False)
    assert verify.expected_outputs(base_config(only_MANE_select=False), '/refs')[0] == \
        ('data/transcripts_reference_113.csv', False)


def test_report(tmp_path, capsys):
    a, b = plain(tmp_path / 'a', TEXT), plain(tmp_path / 'b', TEXT.replace(b'3', b'4'))
    results = [compare.compare('same.tsv', a, a), compare.compare('diff.tsv', a, b),
               compare.compare('gone.tsv', a, str(tmp_path / 'nope'))]
    times = {'legacy': (12.4, 30.0, 0), 'fast': (1.2, 2.0, 1)}
    verify.write_report(str(tmp_path / 'report.tsv'), results, times)
    lines = (tmp_path / 'report.tsv').read_text().splitlines()
    assert lines[:3] == ['# legacy: 12 s wall, 30 CPU-s, exit 0', '# fast: 1 s wall, 2 CPU-s, exit 1',
                         'file\tlegacy_sha256\tfast_sha256\tresult\tfirst_difference']
    assert lines[4].split('\t') == ['diff.tsv', results[1].sha_a, results[1].sha_b, 'DIFFERENT', 'line 3']
    assert lines[5].split('\t') == ['gone.tsv', results[2].sha_a, '-', 'MISSING', 'only in the legacy arm']
    verify.print_report(results, times)
    out = capsys.readouterr().out
    assert 'diff.tsv, line 3:\n  legacy: line 3\n  fast:   line 4' in out
    assert 'FAILED (exit 1)' in out
