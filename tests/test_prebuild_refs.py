"""container/prebuild_refs.py on a tiny synthetic Ensembl release (tests/tiny_ensembl.py).

Runs in the ProHap environment (the reference scripts need gffutils and pandas):

    pytest tests/test_prebuild_refs.py
"""
import gzip
import hashlib
import shutil
import sqlite3
import stat
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tiny_ensembl as tiny  # noqa: E402

pytest.importorskip('gffutils')

ROOT = Path(__file__).resolve().parent.parent
PREBUILD = ROOT / 'container' / 'prebuild_refs.py'
A = tiny.ANNOTATION


def prebuild(downloads, out, *extra):
    return subprocess.run([sys.executable, str(PREBUILD), '--release', tiny.RELEASE, '--downloads', str(downloads),
                           '--out', str(out), '--chromosomes', ','.join(tiny.CHROMOSOMES), '-j', '2', *extra],
                          capture_output=True, text=True)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.fixture(scope='module')
def built(tmp_path_factory):
    base = tmp_path_factory.mktemp('prebuild')
    downloads = Path(tiny.write_downloads(base / 'downloads'))
    proc = prebuild(downloads, base / 'refs')
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return downloads, base / 'refs' / ('ensembl.v' + tiny.RELEASE)


def test_file_set(built):
    _, refs = built
    expected = {A + '.db', A + '_chr21.db', A + '_chr22.db', 'transcripts_reference_113.csv', 'total_cdnas_113.fa',
                'ensembl_reference_proteinDB_113_tagged.fa', 'crap_tagged.fa.gz', 'release.txt', 'SOURCES.sha256',
                'MANIFEST.sha256'}
    assert {p.name for p in refs.iterdir()} == expected
    assert not (refs.parent / 'ensembl.v113.partial').exists()


def test_manifest_and_sources(built):
    downloads, refs = built
    manifest = [line.split('  ', 1) for line in (refs / 'MANIFEST.sha256').read_text().splitlines()]
    assert [name for _, name in manifest] == sorted(p.name for p in refs.iterdir() if p.name != 'MANIFEST.sha256')
    for digest, name in manifest:
        assert sha256(refs / name) == digest, name
    sources = dict(reversed(line.split('  ', 1)) for line in (refs / 'SOURCES.sha256').read_text().splitlines())
    expected = {p.name: sha256(p) for p in downloads.iterdir()}
    expected['crap.fasta'] = sha256(ROOT / 'crap.fasta')
    assert sources == expected


def test_contents(built):
    downloads, refs = built
    assert (refs / 'release.txt').read_text() == '113\n'
    assert (refs / 'transcripts_reference_113.csv').read_text() == tiny.MANE_LIST
    # merge_cdnas_fasta: ncRNA first, then cDNA
    plain = {p.name: gzip.decompress(p.read_bytes()) for p in downloads.iterdir()}
    assert (refs / 'total_cdnas_113.fa').read_bytes() == \
        plain['Homo_sapiens.GRCh38.ncrna.fa.gz'] + plain['Homo_sapiens.GRCh38.cdna.all.fa.gz']
    headers = [line for line in (refs / 'ensembl_reference_proteinDB_113_tagged.fa').read_text().splitlines()
               if line.startswith('>')]
    assert [h.split('|')[:2] for h in headers] == [['>generic_ensref', 'ENST00000000211'],
                                                   ['>generic_ensref', 'ENST00000000212'],
                                                   ['>generic_ensref', 'ENST00000000221'],
                                                   ['>generic_ensref', 'ENST00000000991']]
    with gzip.open(refs / 'crap_tagged.fa.gz', 'rt') as fh:
        tagged = [line for line in fh if line.startswith('>')]
    assert len(tagged) == (ROOT / 'crap.fasta').read_text().count('>')
    assert all(h.startswith('>generic_') for h in tagged)


def test_databases(built):
    _, refs = built
    for name, seqids in [(A + '.db', {'21', '22', 'KI270728.1'}), (A + '_chr21.db', {'21'}),
                         (A + '_chr22.db', {'22'})]:
        con = sqlite3.connect('file:%s?mode=ro' % (refs / name), uri=True)
        try:
            assert {r[0] for r in con.execute('SELECT DISTINCT seqid FROM features')} == seqids
            # read from a read-only directory: anything but WAL
            assert con.execute('PRAGMA journal_mode').fetchone()[0].lower() != 'wal'
        finally:
            con.close()


def test_permissions(built):
    _, refs = built
    assert stat.S_IMODE(refs.stat().st_mode) == 0o755
    for p in refs.iterdir():
        assert stat.S_IMODE(p.stat().st_mode) == 0o644, p.name


def test_refuses_to_overwrite(built):
    downloads, refs = built
    proc = prebuild(downloads, refs.parent)
    assert proc.returncode != 0
    assert 'already exists' in proc.stderr


def test_missing_download(built, tmp_path):
    downloads, _ = built
    partial = tmp_path / 'downloads'
    shutil.copytree(downloads, partial)
    (partial / 'Homo_sapiens.GRCh38.pep.all.fa.gz').unlink()
    proc = prebuild(partial, tmp_path / 'refs')
    assert proc.returncode != 0
    assert 'Homo_sapiens.GRCh38.pep.all.fa.gz' in proc.stderr
    assert not (tmp_path / 'refs' / 'ensembl.v113').exists()


def test_failed_step_leaves_no_release_dir(built, tmp_path):
    downloads, _ = built
    # chromosome 7 has no records, so split_gtf's grep fails, as the upstream rule would
    proc = subprocess.run([sys.executable, str(PREBUILD), '--release', '113', '--downloads', str(downloads),
                           '--out', str(tmp_path / 'refs'), '--chromosomes', '21,7', '-j', '2'],
                          capture_output=True, text=True)
    assert proc.returncode != 0
    assert 'split_gtf chr7 failed' in proc.stderr
    assert not (tmp_path / 'refs' / 'ensembl.v113').exists()
