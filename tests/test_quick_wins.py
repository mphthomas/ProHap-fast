"""The quick wins in src/ (IMPLEMENTATION_NOTES section 1) against legacy/src, function by function, on small
crafted inputs. P1-P3 (get_haplotypes.py, process_haplotypes.py) went with the rewrite of prohap.py in 0.2.0, whose
modules have their own tests (test_sweep, test_haplotypes, test_proteins, test_prohap_fast); vcf_reader.py (P4)
remains for src/provar.py. tests/legacy_vs_fast.sh is the real check; these run in the unit CI, in a few seconds:

    pytest tests/test_quick_wins.py
"""
import importlib
import io
import os
import random
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import pytest

pd = pytest.importorskip('pandas')

ROOT = Path(__file__).resolve().parent.parent
TREES = {'legacy': ROOT / 'legacy' / 'src', 'fast': ROOT / 'src'}
VCF_COLUMNS = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']


@contextmanager
def modules_of(arm):
    """Import `modules.<name>` from one tree. `modules` is a namespace package in both trees, so only one tree may be
    on sys.path at a time and its imports are dropped afterwards."""
    def drop():
        for name in [n for n in sys.modules if n == 'modules' or n.startswith('modules.')]:
            del sys.modules[name]
    drop()
    sys.path.insert(0, str(TREES[arm]))
    try:
        yield lambda name: importlib.import_module('modules.' + name)
    finally:
        sys.path.remove(str(TREES[arm]))
        drop()


# ---------------------------------------------------------------------------------------------------------------
# P4: vcf_reader.py collects each transcript's VCF lines in a list and joins them once


class FakeDB:
    """The one gffutils call parse_vcf makes: the exons of a transcript, ordered by start."""
    def __init__(self, exons):
        self.exons = exons

    def children(self, transcript, featuretype, order_by):
        assert (featuretype, order_by) == ('exon', 'start')
        return [SimpleNamespace(start=s, end=e) for s, e in self.exons[transcript.id]]


def test_parse_vcf_matches_legacy(tmp_path):
    rng = random.Random(3)
    exons = {'T1': [(1000, 1100), (1500, 1600)], 'T2': [(1050, 1080), (1550, 1700)], 'T3': [(1520, 1530)],
             'T4': [(5000, 5100)], 'T5': [(9000, 9100)]}
    transcripts = [SimpleNamespace(id=t, start=e[0][0], end=e[-1][1]) for t, e in exons.items()]
    samples = ['S%d' % i for i in range(1, 9)]
    lines = ['##fileformat=VCFv4.2', '##contig=<ID=22>', '#' + '\t'.join(VCF_COLUMNS + samples)]
    positions = [500, 1000, 1001, 1060, 1098, 1100, 1300, 1495, 1525, 1525, 1560, 1650, 3000, 5050, 5100, 7000]
    for n, pos in enumerate(positions):
        ref, alt = rng.choice([('A', 'G'), ('C', 'T'), ('ACGTA', 'A'), ('G', 'GTT'), ('T', 'C,G'), ('N', 'A')])
        vcf_id = rng.choice(['.', '.', 'rs%d' % n, '2e5', 'NA'])
        af = rng.choice(['0.01', '0.2', '0.5'])
        gts = ['%d|%d' % (rng.random() < 0.3, rng.random() < 0.3) for _ in samples]
        lines.append('\t'.join(['22', str(pos), vcf_id, ref, alt, '.', 'PASS', 'AC=1;AN=16;AF=' + af, 'GT'] + gts))
    text = '\n'.join(lines) + '\n'

    written = {}
    for arm in ('legacy', 'fast'):
        out = tmp_path / arm
        out.mkdir()
        with modules_of(arm) as module:
            colnames = module('vcf_reader').parse_vcf(transcripts, io.StringIO(text), FakeDB(exons), 0.05, str(out))
        written[arm] = colnames, {p.name: p.read_bytes() for p in sorted(out.iterdir())}
    assert written['fast'] == written['legacy']
    tables = written['legacy'][1]
    assert sorted(tables) == ['T%d.tsv' % i for i in range(1, 6)]
    assert sum(len(t.splitlines()) - 1 for t in tables.values()) > len(exons)    # lines shared by transcripts


# ---------------------------------------------------------------------------------------------------------------
# D3: merge_duplicate_seq.py merges by sequence hash in a dict and sorts once


def test_merge_duplicate_seq_matches_legacy(tmp_path):
    """Both scripts on 60 records of 15 sequences, with every kind of header the script parses. The record order is
    by hash(), so the scripts run with PYTHONHASHSEED=0 as in the workflow."""
    rng = random.Random(5)
    sequences = [''.join(rng.choice('ACDEFGHIKLMN') for _ in range(rng.randint(5, 30))) for _ in range(15)]
    headers = [
        'generic_ensref|ENSP{i:05d}|matching_proteins:ENST{i:05d} position_within_protein:{r}',
        'generic_enshap|enshap_22_{i:x}|matching_proteins:haplo_chr22_{i:x};haplo_chr22_{j:x} start:{r} '
        'reading_frame:{r};{s}',
        'generic_var|var_{i}|matching_proteins:var_{i} start:0 reading_frame:{r} position_within_protein:{s}',
        'generic_cont|CONT{i}|contaminant {i}',
        'sp|P{i:05d}|Swiss-Prot entry {i}',
    ]
    fasta = tmp_path / 'input.fa'
    with open(fasta, 'w') as out:
        for i in range(60):
            header = headers[i % 5].format(i=i, j=i + 100, r=rng.randint(0, 2), s=rng.randint(0, 2))
            out.write('>' + header + '\n' + rng.choice(sequences) + '\n')

    env = dict(os.environ, PYTHONHASHSEED='0', PYTHONNOUSERSITE='1')
    outputs = {}
    for arm in ('legacy', 'fast'):
        result = tmp_path / (arm + '.fa')
        subprocess.run([sys.executable, str(TREES[arm] / 'merge_duplicate_seq.py'), '-i', str(fasta), '-o',
                        str(result)], env=env, check=True)
        outputs[arm] = result.read_text()
    assert outputs['fast'] == outputs['legacy']
    assert outputs['legacy'].count('>') == len(set(sequences))
