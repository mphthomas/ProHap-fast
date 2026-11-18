"""How the Snakefile finds its references: refs_dir, the release checks, the upstream keys, the transcript list and
the contaminants. Dry runs against a tiny synthetic release (tests/tiny_ensembl.py) built by prebuild_refs.py.

Runs in the ProHap environment (snakemake 7, gffutils):

    pytest tests/test_snakefile_refs.py
"""
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tiny_ensembl as tiny  # noqa: E402

pytest.importorskip('gffutils')

ROOT = Path(__file__).resolve().parent.parent
SNAKEMAKE = Path(sys.executable).parent / 'snakemake'
if not SNAKEMAKE.exists():
    pytest.skip('snakemake is not installed next to ' + sys.executable, allow_module_level=True)
REMOVED_RULES = {'download_vcf', 'download_gtf', 'download_cdnas_fasta', 'download_reference_proteome',
                 'merge_cdnas_fasta', 'parse_gtf_whole', 'split_gtf', 'parse_gtf_chromosome',
                 'reference_filter_format'}


@pytest.fixture(scope='module')
def refs_root(tmp_path_factory):
    base = tmp_path_factory.mktemp('refs')
    downloads = tiny.write_downloads(base / 'downloads')
    proc = subprocess.run([sys.executable, str(ROOT / 'container' / 'prebuild_refs.py'), '--release', tiny.RELEASE,
                           '--downloads', str(downloads), '--out', str(base / 'refs'),
                           '--chromosomes', ','.join(tiny.CHROMOSOMES), '-j', '2'], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return base / 'refs'


@pytest.fixture
def workdir(tmp_path):
    (tmp_path / 'samples.tsv').write_text('Sample name\tSex\tPopulation code\tSuperpopulation code\n'
                                          'S1\tfemale\tFIN\tEUR\n')
    return tmp_path


def base_config(workdir, refs_root):
    return {
        'implementation': 'fast', 'refs_dir': str(refs_root), 'ensembl_release': 113, 'chromosomes': ['21', '22'],
        'contaminants_fasta': str(ROOT / 'crap.fasta'), 'included_transcript_biotypes': '', 'only_MANE_select': True,
        'final_fasta_file': 'results/db.fa', 'fasta_simplify_headers': True, 'use_ProHap': True,
        'phased_local_path': str(workdir) + '/', 'phased_vcf_file_name': 'chr{chr}.vcf.gz',
        'sample_metadata_file': str(workdir / 'samples.tsv'), 'working_dir_name_haplo': 'haplotypes_tmp',
        'phased_min_af': 0.01, 'phased_af_field': 'AF', 'haplo_min_freq': 0, 'haplo_min_count': 0,
        'x_par1_to': 2781479, 'x_par2_from': 155701383, 'haplo_require_start': 1, 'haplo_ignore_UTR': 1,
        'haplo_skip_start_lost': 1, 'max_cores': 1, 'haplo_fasta_file': 'results/haplo.fa.gz',
        'haplo_cdna_file': '', 'haplo_table_file': 'results/haplo.tsv.gz', 'use_ProVar': False, 'variant_vcf': {},
        'working_dir_name_var': 'variants_tmp', 'var_require_start': 1, 'var_fasta_file': '', 'var_cdna_file': '',
        'var_table_file': '', 'add_existing_haplo': False, 'haplo_added_table': '', 'haplo_added_fasta': '',
    }


class DryRun:
    def __init__(self, proc):
        self.returncode = proc.returncode
        self.output = proc.stdout + proc.stderr
        self.rules = {line.split()[1].rstrip(':') for line in self.output.splitlines()
                      if line.startswith(('rule ', 'localrule '))}

    def command(self, rule):
        """The shell command printed for the first job of `rule`."""
        lines = self.output.splitlines()
        start = next(i for i, line in enumerate(lines) if line.split()[:2] in (['rule', rule + ':'],
                                                                                ['localrule', rule + ':']))
        for line in lines[start + 1:]:
            if line.startswith(('rule ', 'localrule ', 'Job stats')):
                break
            if 'python3' in line or line.startswith('cat '):
                return line
        raise AssertionError('no command for ' + rule)


def dry_run(workdir, config, **env):
    (workdir / 'config.yaml').write_text(yaml.safe_dump(config))
    environ = {k: v for k, v in os.environ.items() if k not in ('ENSEMBL_RELEASE', 'PROHAP_REFS')}
    environ.update(env)
    proc = subprocess.run([str(SNAKEMAKE), '--snakefile', str(ROOT / 'Snakefile'), '--directory', str(workdir),
                           '--cores', '1', '-n', '-p'], capture_output=True, text=True, env=environ, cwd=workdir)
    return DryRun(proc)


def ref(refs_root, name):
    return str(refs_root / 'ensembl.v113' / name)


def test_defaults_read_everything_from_refs(workdir, refs_root):
    run = dry_run(workdir, base_config(workdir, refs_root))
    assert run.returncode == 0, run.output
    assert not run.rules & (REMOVED_RULES | {'default_transcript_list', 'contaminants_fix_headers'})
    cmd = run.command('compute_haplotypes')
    assert '-db "%s"' % ref(refs_root, tiny.ANNOTATION + '_chr21.db') in cmd or \
        '-db "%s"' % ref(refs_root, tiny.ANNOTATION + '_chr22.db') in cmd
    assert '-transcripts "%s"' % ref(refs_root, 'transcripts_reference_113.csv') in cmd
    assert '-cdna "%s"' % ref(refs_root, 'total_cdnas_113.fa') in cmd
    cmd = run.command('reference_remove_stop')
    assert '-i %s' % ref(refs_root, 'ensembl_reference_proteinDB_113_tagged.fa') in cmd
    assert '-tr %s' % ref(refs_root, 'transcripts_reference_113.csv') in cmd
    assert '-db %s' % ref(refs_root, tiny.ANNOTATION + '.db') in run.command('simpify_fasta_headers')
    assert 'warning' not in run.output.lower()


@pytest.mark.parametrize('listed', ['data/transcripts_reference_113.csv', './data/transcripts_reference_113.csv'])
def test_upstream_reference_list_path_means_the_prebuilt_list(workdir, refs_root, listed):
    config = dict(base_config(workdir, refs_root), custom_transcript_list=listed)
    run = dry_run(workdir, config)
    assert run.returncode == 0, run.output
    assert '-transcripts "%s"' % ref(refs_root, 'transcripts_reference_113.csv') in run.command('compute_haplotypes')


def test_empty_list_means_every_transcript_of_the_biotypes(workdir, refs_root):
    run = dry_run(workdir, dict(base_config(workdir, refs_root), custom_transcript_list=''))
    assert run.returncode == 0, run.output
    assert 'get_transcript_list' in run.rules
    assert '-i %s' % ref(refs_root, tiny.ANNOTATION + '.db') in run.command('get_transcript_list')
    assert '-transcripts "data/included_transcripts.csv"' in run.command('compute_haplotypes')


def test_own_list_is_used_as_is(workdir, refs_root):
    (workdir / 'mine.csv').write_text('chromosome,transcriptID\n21,ENST00000000211\n')
    run = dry_run(workdir, dict(base_config(workdir, refs_root), custom_transcript_list='mine.csv'))
    assert run.returncode == 0, run.output
    assert '-transcripts "mine.csv"' in run.command('compute_haplotypes')


def test_without_mane_the_reference_list_is_built_in_the_work_dir(workdir, refs_root):
    run = dry_run(workdir, dict(base_config(workdir, refs_root), only_MANE_select=False))
    assert run.returncode == 0, run.output
    cmd = run.command('default_transcript_list')
    assert '-i %s' % ref(refs_root, 'ensembl_reference_proteinDB_113_tagged.fa') in cmd
    assert '-annot %s' % ref(refs_root, tiny.ANNOTATION + '.db') in cmd
    assert '-MANE 0 -o data/transcripts_reference_113.csv' in cmd
    assert '-tr' not in run.command('reference_remove_stop')
    assert '-transcripts "data/transcripts_reference_113.csv"' in run.command('compute_haplotypes')


def test_contaminants(workdir, refs_root):
    shutil.copy(ROOT / 'crap.fasta', workdir / 'same.fasta')
    run = dry_run(workdir, dict(base_config(workdir, refs_root), contaminants_fasta='same.fasta'))
    assert run.returncode == 0, run.output
    assert 'contaminants_fix_headers' not in run.rules

    (workdir / 'other.fasta').write_text('>sp|P00000|TEST_HUMAN test\nMKKK\n')
    run = dry_run(workdir, dict(base_config(workdir, refs_root), contaminants_fasta='other.fasta'))
    assert run.returncode == 0, run.output
    assert 'contaminants_fix_headers' in run.rules
    assert run.command('contaminants_fix_headers').endswith('-i other.fasta -o data/fasta/crap_tagged.fa.gz -t _cont')


def test_refs_dir_defaults_to_prohap_refs(workdir, refs_root):
    config = base_config(workdir, refs_root)
    del config['refs_dir']
    run = dry_run(workdir, config, PROHAP_REFS=str(refs_root))
    assert run.returncode == 0, run.output
    assert ref(refs_root, 'total_cdnas_113.fa') in run.command('compute_haplotypes')


@pytest.mark.parametrize('key, value', [
    ('Ensembl_FTP_URL', 'ftp.ensembl.org/pub/release-113/'),
    ('annotationFilename', tiny.ANNOTATION),
    ('phased_FTP_URL', 'http://example.org/vcf/'),
])
def test_consistent_upstream_keys_warn(workdir, refs_root, key, value):
    run = dry_run(workdir, dict(base_config(workdir, refs_root), **{key: value}))
    assert run.returncode == 0, run.output
    assert "ProHap-fast config warning: '%s' is ignored" % key in run.output


@pytest.mark.parametrize('key, value', [('Ensembl_FTP_URL', ''), ('annotationFilename', ''), ('phased_FTP_URL', '')])
def test_empty_upstream_keys_are_silent(workdir, refs_root, key, value):
    run = dry_run(workdir, dict(base_config(workdir, refs_root), **{key: value}))
    assert run.returncode == 0, run.output
    assert 'warning' not in run.output.lower()


@pytest.mark.parametrize('changes, message', [
    ({'Ensembl_FTP_URL': 'ftp.ensembl.org/pub/release-114/'}, "does not name GRCh38 release 113"),
    ({'Ensembl_FTP_URL': 'ftp.ensembl.org/pub/current/'}, "does not name GRCh38 release 113"),
    ({'Ensembl_FTP_URL': 'ftp.ensembl.org/pub/grch37/release-113/'}, "does not name GRCh38 release 113"),
    ({'annotationFilename': 'Homo_sapiens.GRCh38.114.chr_patch_hapl_scaff'}, "delete 'annotationFilename'"),
    ({'ensembl_release': 114}, "no references for Ensembl release 114 (config 'ensembl_release')"),
    ({'ensembl_release': 'latest'}, "must be a release number"),
    ({'chromosomes': ['21', 'X']}, "needs it"),
    ({'contaminants_fasta': 'missing.fasta'}, "'contaminants_fasta' is 'missing.fasta', which is not a file"),
])
def test_errors(workdir, refs_root, changes, message):
    run = dry_run(workdir, dict(base_config(workdir, refs_root), **changes))
    assert run.returncode != 0
    assert message in run.output, run.output
    assert 'ENSEMBL_RELEASE is not set, so this is not running in a prohap-fast image.' in run.output


def test_missing_chromosome_names_the_key(workdir, refs_root):
    run = dry_run(workdir, dict(base_config(workdir, refs_root), chromosomes=['21', 'X']))
    assert "config 'chromosomes' (chromosome X) needs it" in run.output
    assert ref(refs_root, tiny.ANNOTATION + '_chrX.db') in run.output


def test_wrong_release_names_the_image_release(workdir, refs_root):
    run = dry_run(workdir, dict(base_config(workdir, refs_root), ensembl_release=114), ENSEMBL_RELEASE='113')
    assert run.returncode != 0
    assert 'has release(s) 113' in run.output
    assert 'prohap-fast:<version>-ensembl.v114' in run.output
    assert 'This image carries the references of Ensembl release 113.' in run.output


def test_release_txt_must_match(workdir, refs_root, tmp_path_factory):
    other = tmp_path_factory.mktemp('relabelled')
    shutil.copytree(refs_root / 'ensembl.v113', other / 'ensembl.v113')
    release_txt = other / 'ensembl.v113' / 'release.txt'
    release_txt.chmod(release_txt.stat().st_mode | stat.S_IWUSR)
    release_txt.write_text('114\n')
    run = dry_run(workdir, dict(base_config(workdir, refs_root), refs_dir=str(other)))
    assert run.returncode != 0
    assert "release.txt says release '114', but config 'ensembl_release' is 113" in run.output


def test_refs_dir_with_a_space(workdir, refs_root, tmp_path_factory):
    spaced = tmp_path_factory.mktemp('with space')
    run = dry_run(workdir, dict(base_config(workdir, refs_root), refs_dir=str(spaced)))
    assert run.returncode != 0
    assert 'contains a space or a shell character' in run.output
