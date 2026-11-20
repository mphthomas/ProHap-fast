"""D4: merge_duplicate_seq.py and remove_UTR_seq.py read the FASTA record by record (src/modules/fasta_stream.py).

modules.fasta_stream.records against legacy read_fasta on odd FASTA files, then both scripts against the reference
scripts, as scripts, with PYTHONHASHSEED=0: records spread over the hash partitions, every kind of header, UTR
fragments, a repeated accession (merge: the dict way; remove_UTR: the reference) and an empty file. The dupseq
fixture of tests/legacy_vs_fast.sh (medium tier) runs both at 400,000 records.

    pytest tests/test_fasta_stream.py
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
sys.path.insert(0, str(ROOT / 'src'))
from modules import fasta_stream  # noqa: E402
from modules.legacy_import import legacy  # noqa: E402

ENV = dict(os.environ, PYTHONHASHSEED='0', PYTHONNOUSERSITE='1')

HEADERS = [
    'generic_ensref|ENSP{i:05d}|matching_proteins:ENST{i:05d} position_within_protein:{r}',
    'generic_ensref|ENSP{i:05d}.{r}|matching_proteins:ENST{i:05d}',
    'generic_enshap|enshap_22_{i:x}|matching_proteins:haplo_chr22_{i:x};haplo_chr22_{j:x} start:{r} '
    'reading_frame:{r};{s} position_within_protein:{s}',
    'generic_enshap|enshap_22_{i:x}_3UTR_1|matching_proteins:haplo_chr22_{i:x} start:0 reading_frame:{r}',
    'generic_var|var_{i}|matching_proteins:var_{i} start:0 reading_frame:{r} position_within_protein:{s}',
    'generic_cont|CONT{i}|contaminant {i}',
    'sp|P{i:05d}|Swiss-Prot entry {i}',
    'generic_cont|NODESC{i}',
    'generic_cont|ACC{i} merged description {i}',
    'PLAIN{i} no bars here',
    'BARE{i}',
]


def fasta_text(rng, n, n_seqs, wrap=False, repeat=None):
    sequences = [''.join(rng.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(rng.randint(3, 40))) for _ in range(n_seqs)]
    out = []
    for i in range(n):
        header = HEADERS[i % len(HEADERS)].format(i=i, j=i + 1000, r=rng.randint(0, 2), s=rng.randint(0, 2))
        if repeat is not None and i == repeat[1]:
            header = HEADERS[repeat[0] % len(HEADERS)].format(i=repeat[0], j=repeat[0] + 1000, r=0, s=0)
        seq = rng.choice(sequences)
        body = '\n'.join(seq[k:k + 10] for k in range(0, len(seq), 10)) if wrap else seq
        out.append('>' + header + '\n' + body + '\n')
    return ''.join(out)


@pytest.mark.parametrize('variant', ['plain', 'wrapped', 'no_final_newline', 'not_starting_with_header', 'empty'])
def test_records_as_read_fasta(tmp_path, variant):
    rng = random.Random(1)
    text = fasta_text(rng, 50, 20, wrap=(variant == 'wrapped'))
    if variant == 'no_final_newline':
        text = text.rstrip('\n')
    elif variant == 'not_starting_with_header':
        text = 'stray line\n' + text
    elif variant == 'empty':
        text = ''
    path = tmp_path / 'in.fa.gz'
    with gzip.open(path, 'wt') as fh:
        fh.write(text)
    for truncate in (False, True):
        expected = legacy.common.read_fasta(str(path), truncate)
        streamed = list(fasta_stream.records(str(path), truncate))
        if len({k for k, _ in streamed}) == len(streamed):
            assert dict(streamed) == expected and [k for k, _ in streamed] == list(expected)
        else:       # truncated accessions repeat (ENSP00001.0 and ENSP00001.1): read_fasta keeps the last value
            assert dict(streamed) == expected


def test_accessions_repeated():
    a = fasta_stream.Accessions()
    for x in ('a', 'b', 'c'):
        a.add(x)
    assert not a.repeated() and len(a) == 3
    a.add('b')
    assert a.repeated()
    assert not fasta_stream.Accessions().repeated()


def run(arm, script, src, dst):
    return subprocess.run([sys.executable, str(TREES[arm] / script), '-i', str(src), '-o', str(dst)],
                          capture_output=True, text=True, env=ENV)


def read(path):
    data = Path(path).read_bytes()
    return gzip.decompress(data) if data[:2] == b'\x1f\x8b' else data


def both(tmp_path, text, gz=True):
    src = tmp_path / ('in.fa.gz' if gz else 'in.fa')
    if gz:
        with gzip.open(src, 'wt') as fh:
            fh.write(text)
    else:
        src.write_text(text)
    out = {}
    for arm in ('legacy', 'fast'):
        merged = tmp_path / ('%s_nodupl.fa.gz' % arm)
        final = tmp_path / ('%s_final.fa' % arm)
        p1 = run(arm, 'merge_duplicate_seq.py', src, merged)
        assert p1.returncode == 0, p1.stderr
        p2 = run(arm, 'remove_UTR_seq.py', merged, final)
        assert p2.returncode == 0, p2.stderr
        out[arm] = (read(merged), read(final), p1.stdout + p1.stderr, p2.stderr)
    assert out['fast'][:2] == out['legacy'][:2]
    assert not list(tmp_path.glob('merge_duplicate_seq.*')), 'temporary directory left behind'
    return out


@pytest.mark.parametrize('seed,n,n_seqs,gz', [(1, 60, 15, True), (2, 3000, 900, True), (3, 2000, 2000, False)])
def test_scripts_match_the_reference(tmp_path, seed, n, n_seqs, gz):
    out = both(tmp_path, fasta_text(random.Random(seed), n, n_seqs), gz)
    merged = out['legacy'][0].decode()
    assert merged.count('>') >= min(n_seqs, 100) // 2
    assert 'prot_0|' in merged and 'UTR' in merged
    assert 'merging in memory' not in out['fast'][2] and FALLBACK not in out['fast'][3]


def test_repeated_accession(tmp_path):
    """merge_duplicate_seq.py merges the dict way; remove_UTR_seq.py gets numbered records and streams."""
    out = both(tmp_path, fasta_text(random.Random(4), 200, 50, repeat=(5, 120)))
    assert 'merging in memory' in out['fast'][2]
    assert FALLBACK not in out['fast'][3]


def test_remove_utr_repeated_accession_falls_back(tmp_path):
    rng = random.Random(5)
    rec = lambda acc, seq: ('>generic_enshap|%s|position_within_protein:0 start:0 matching_proteins:h1 '
                            'reading_frame:0 split_sequences:x\n%s\n' % (acc, seq))
    src = tmp_path / 'in.fa'
    src.write_text(rec('prot_0', 'MAAA') + rec('prot_1', 'MCCC') + rec('prot_0', 'MDDD'))
    outs = {}
    for arm in ('legacy', 'fast'):
        proc = run(arm, 'remove_UTR_seq.py', src, tmp_path / (arm + '.fa'))
        assert proc.returncode == 0, proc.stderr
        outs[arm] = (tmp_path / (arm + '.fa')).read_text()
        if arm == 'fast':
            assert FALLBACK in proc.stderr
    assert outs['fast'] == outs['legacy'] and outs['legacy'].count('>') == 2


def test_empty_input(tmp_path):
    out = both(tmp_path, '')
    assert out['legacy'][:2] == (b'', b'')
