import hashlib
import os
import re
import shlex
import sys

from snakemake.exceptions import WorkflowError

configfile: "config.yaml"

# Config `implementation` picks the script directory. Every rule, path and reference is shared, so the two arms
# differ only in which scripts run. legacy/src is the frozen reference (upstream cafcc51 + legacy/legacy-fixes.patch)
# and src/ the fast implementation. Absolute paths let the workflow run with --directory pointing anywhere.
LEGACY_SRC = os.path.join(workflow.basedir, 'legacy', 'src')
FAST_SRC = os.path.join(workflow.basedir, 'src')
IMPLEMENTATION = config.get('implementation', 'fast')
if IMPLEMENTATION not in ('fast', 'legacy'):
    raise ValueError("config 'implementation' must be 'fast' or 'legacy', not " + repr(IMPLEMENTATION))
SRC = FAST_SRC if IMPLEMENTATION == 'fast' else LEGACY_SRC

# macOS starts multiprocessing workers with 'spawn', which breaks the closures ProHap hands to Pool. There, every
# script runs through a wrapper that selects 'fork', the default on Linux.
if sys.platform == 'darwin':
    PYTHON = 'python3 ' + shlex.quote(os.path.join(workflow.basedir, 'tools', 'run_fork.py'))
else:
    PYTHON = 'python3'

def script(name, src=None):
    """Shell command that runs one of ProHap's scripts from SRC, or from src (LEGACY_SRC for the ProVar rules)."""
    return PYTHON + ' ' + shlex.quote(os.path.join(src or SRC, name))

# shell.prefix() replaces Snakemake's default prefix, bash strict mode, so it is repeated here. The fixed hash seed
# keeps FASTA record order reproducible, and the user's site-packages must not shadow the environment's packages.
shell.prefix("set -euo pipefail; export PYTHONHASHSEED=0 PYTHONNOUSERSITE=1; ")

VARIANT_VCF_FILES = config['variant_vcf']
CHROMOSOMES = config['chromosomes']

# ------------------------------------ references ------------------------------------
# Everything derived from the Ensembl release is pre-built once by container/prebuild_refs.py into
# <refs_dir>/ensembl.v<release>/ and read from there by both arms: nothing is downloaded and nothing is copied into
# the work directory. ensembl_release is the only release setting; the annotation name and the reference transcript
# list follow from it.
ENSEMBL_RELEASE = str(config['ensembl_release'])
ANNOTATION = 'Homo_sapiens.GRCh38.' + ENSEMBL_RELEASE + '.chr_patch_hapl_scaff'
REFS_ROOT = os.path.abspath(os.path.expanduser(str(config.get('refs_dir') or os.environ.get('PROHAP_REFS') or '/opt/prohap/refs')))
REFS = os.path.join(REFS_ROOT, 'ensembl.v' + ENSEMBL_RELEASE)
# Set in the release image; names the release whose references the image carries.
IMAGE_RELEASE = os.environ.get('ENSEMBL_RELEASE', '')

def config_error(message):
    where = ('This image carries the references of Ensembl release ' + IMAGE_RELEASE + '.') if IMAGE_RELEASE \
        else 'ENSEMBL_RELEASE is not set, so this is not running in a prohap-fast image.'
    raise WorkflowError(message + '\n' + where)

def config_warning(message):
    print('ProHap-fast config warning: ' + message, file=sys.stderr)

if not re.fullmatch(r'[0-9]+', ENSEMBL_RELEASE):
    config_error("config 'ensembl_release' must be a release number such as 113, not " + repr(config['ensembl_release']))
if re.search(r'[^\w./+@%:,=-]', REFS):
    config_error("the references directory " + repr(REFS) + " (config 'refs_dir', else $PROHAP_REFS) contains a space or "
                 "a shell character; use a path without them")

# Upstream release and download settings. References now come only from REFS, so these do nothing: accepted with a
# warning while they agree with ensembl_release, an error when they name something else.
_ftp_url = str(config.get('Ensembl_FTP_URL') or '')
if _ftp_url:
    if re.findall(r'release-([0-9]+)', _ftp_url) != [ENSEMBL_RELEASE] or 'grch37' in _ftp_url.lower():
        config_error("config 'Ensembl_FTP_URL' is " + repr(_ftp_url) + ", which does not name GRCh38 release "
                     + ENSEMBL_RELEASE + " of config 'ensembl_release'. ProHap-fast downloads nothing: delete "
                     "'Ensembl_FTP_URL' and set 'ensembl_release' to the release you want.")
    config_warning("'Ensembl_FTP_URL' is ignored; the release " + ENSEMBL_RELEASE + " references are read from " + REFS
                   + ". Delete the line.")
_annotation = str(config.get('annotationFilename') or '')
if _annotation:
    if _annotation != ANNOTATION:
        config_error("config 'annotationFilename' is " + repr(_annotation) + ", but config 'ensembl_release' "
                     + ENSEMBL_RELEASE + " means " + repr(ANNOTATION) + ". The annotation follows from "
                     "'ensembl_release': delete 'annotationFilename'.")
    config_warning("'annotationFilename' is ignored; it follows from 'ensembl_release'. Delete the line.")
if str(config.get('phased_FTP_URL') or ''):
    config_warning("'phased_FTP_URL' is ignored; nothing is downloaded, and the phased VCFs are read from "
                   "'phased_local_path'. Delete the line.")

def ref(name):
    return os.path.join(REFS, name)

REF_WHOLE_DB = ref(ANNOTATION + '.db')
REF_CHR_DB = ref(ANNOTATION + '_chr{chr}.db')
REF_CDNA = ref('total_cdnas_' + ENSEMBL_RELEASE + '.fa')
REF_PROTEOME = ref('ensembl_reference_proteinDB_' + ENSEMBL_RELEASE + '_tagged.fa')
REF_MANE_LIST = ref('transcripts_reference_' + ENSEMBL_RELEASE + '.csv')
REF_CONTAMINANTS = ref('crap_tagged.fa.gz')

if not os.path.isdir(REFS):
    _found = sorted(d[len('ensembl.v'):] for d in (os.listdir(REFS_ROOT) if os.path.isdir(REFS_ROOT) else [])
                    if d.startswith('ensembl.v') and os.path.isdir(os.path.join(REFS_ROOT, d)))
    config_error("no references for Ensembl release " + ENSEMBL_RELEASE + " (config 'ensembl_release'): " + REFS
                 + " does not exist. The references directory " + REFS_ROOT + " (config 'refs_dir', else $PROHAP_REFS) "
                 + ("has release(s) " + ', '.join(_found) if _found else "has none")
                 + ". Use the image for release " + ENSEMBL_RELEASE + " (prohap-fast:<version>-ensembl.v"
                 + ENSEMBL_RELEASE + "), or set 'ensembl_release' to the release of your references.")
try:
    with open(ref('release.txt')) as _fh:
        _refs_release = _fh.read().strip()
except OSError:
    _refs_release = None
if _refs_release != ENSEMBL_RELEASE:
    config_error(ref('release.txt') + (" is missing" if _refs_release is None else " says release " + repr(_refs_release))
                 + ", but config 'ensembl_release' is " + ENSEMBL_RELEASE + ". Rebuild the directory with "
                 "container/prebuild_refs.py --release " + ENSEMBL_RELEASE + ".")

# The tagged contaminants in REFS were built from the contaminants FASTA recorded in SOURCES.sha256 (the crap.fasta
# shipped with ProHap). Any other contaminants_fasta is tagged in the work directory, as upstream did.
def _sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()

def _refs_contaminants_sha256():
    try:
        with open(ref('SOURCES.sha256')) as fh:
            entries = [line.rstrip('\n').split('  ', 1) for line in fh]
    except OSError:
        return None
    # every other source is an Ensembl download, named Homo_sapiens.GRCh38.*
    others = [e[0] for e in entries if len(e) == 2 and not e[1].startswith('Homo_sapiens.GRCh38.')]
    return others[0] if len(others) == 1 else None

CONTAMINANTS_FASTA = config['contaminants_fasta']
if not os.path.isfile(CONTAMINANTS_FASTA):
    config_error("config 'contaminants_fasta' is " + repr(CONTAMINANTS_FASTA) + ", which is not a file")
TAGGED_CONTAMINANTS = REF_CONTAMINANTS if _sha256(CONTAMINANTS_FASTA) == _refs_contaminants_sha256() \
    else 'data/fasta/crap_tagged.fa.gz'

# Transcript list. Upstream set custom_transcript_list to data/transcripts_reference_<release>.csv, which the
# default_transcript_list rule built in the work directory; with only_MANE_select that list is in REFS. An absent
# custom_transcript_list means the same list; an empty one means every transcript of included_transcript_biotypes,
# as upstream.
WORKDIR_REFERENCE_LIST = 'data/transcripts_reference_' + ENSEMBL_RELEASE + '.csv'
REFERENCE_LIST = REF_MANE_LIST if config['only_MANE_select'] else WORKDIR_REFERENCE_LIST
if config.get('custom_transcript_list') is None:
    TRANSCRIPT_LIST = REFERENCE_LIST
elif len(config['custom_transcript_list']) == 0:
    TRANSCRIPT_LIST = 'data/included_transcripts.csv'
elif os.path.normpath(config['custom_transcript_list']) == WORKDIR_REFERENCE_LIST:
    TRANSCRIPT_LIST = REFERENCE_LIST
else:
    TRANSCRIPT_LIST = config['custom_transcript_list']

# Fail at parse time, naming the config key, rather than with Snakemake's MissingInputException.
_required = [(REF_WHOLE_DB, "'ensembl_release'"), (REF_PROTEOME, "'ensembl_release'")]
if config['only_MANE_select']:
    _required.append((REF_MANE_LIST, "'only_MANE_select'"))
if TAGGED_CONTAMINANTS == REF_CONTAMINANTS:
    _required.append((REF_CONTAMINANTS, "'contaminants_fasta'"))
if config['use_ProHap'] or config['use_ProVar']:
    _required.append((REF_CDNA, "'ensembl_release'"))
    _required += [(REF_CHR_DB.format(chr=c), "'chromosomes' (chromosome " + str(c) + ")") for c in CHROMOSOMES]
for _path, _key in _required:
    if not os.path.isfile(_path):
        config_error("reference " + _path + " is missing; config " + _key + " needs it. The references for release "
                     + ENSEMBL_RELEASE + " are incomplete: rebuild them with container/prebuild_refs.py.")

WORKING_DIR_NAME_HAPLO = config['working_dir_name_haplo']
WORKING_DIR_NAME_VAR = config['working_dir_name_var']

rule all:
    input:
        final_fasta=expand('{proxy}', proxy=[config['final_fasta_file']] if not config['fasta_simplify_headers'] else ['.'.join(config['final_fasta_file'].split('.')[:-1]) + '_simplified.fasta']),
        var_table=expand('{proxy}', proxy=[config['var_table_file']] if config["use_ProVar"] else []),
        haplo_table=expand('{proxy}', proxy=[config['haplo_table_file']] if config["use_ProHap"] else []),
        var_fasta=expand('{proxy}', proxy=[config['var_fasta_file']] if config["use_ProVar"] else []),
        haplo_fasta=expand('{proxy}', proxy=[config['haplo_fasta_file']] if config["use_ProHap"] else []),

rule get_transcript_list:
    input:
        REF_WHOLE_DB
    output:
        temp("data/included_transcripts.csv")
    params:
        biotypes=config['included_transcript_biotypes']
    shell:
        script('get_transcript_list.py') + " -i {input} -bio {params.biotypes} -o {output}"

# With only_MANE_select the reference transcript list is pre-built in REFS; otherwise it is built here.
if not config['only_MANE_select']:
    rule default_transcript_list:
        input:
            ref_fasta=REF_PROTEOME,
            annot=REF_WHOLE_DB
        output:
            # not temp(): reused across runs, and often given directly as custom_transcript_list
            WORKDIR_REFERENCE_LIST
        params:
            MANE=int(config['only_MANE_select'])
        conda: "envs/prohap.yaml"
        shell:
            script('get_reference_ENST.py') + " -i {input.ref_fasta} -annot {input.annot} -MANE {params.MANE} -o {output}"

rule reference_remove_stop:
    input:
        fasta=REF_PROTEOME,
        tr=expand('{proxy}', proxy=[REF_MANE_LIST] if config['only_MANE_select'] else [])
    output:
        temp("data/fasta/ensembl_reference_proteinDB_" + ENSEMBL_RELEASE + "_clean.fa.gz")
    params:
        tr_filter=('-tr ' + REF_MANE_LIST) if config['only_MANE_select'] else ''
    conda: "envs/prohap.yaml"
    shell:
        script('remove_stop_codons.py') + " -i {input.fasta} -o {output} -min_len 6 {params.tr_filter} "

# The pre-built tagged contaminants come from ProHap's crap.fasta; any other contaminants_fasta is tagged here.
if TAGGED_CONTAMINANTS != REF_CONTAMINANTS:
    rule contaminants_fix_headers:
        input:
            CONTAMINANTS_FASTA
        output:
            TAGGED_CONTAMINANTS
        conda: "envs/prohap.yaml"
        shell:
            script('fasta_format_headers.py') + " -i {input} -o {output} -t _cont"

# ------------------------------------ ProVar rules ------------------------------------

rule split_variant_vcf:
    input:
        lambda wildcards: VARIANT_VCF_FILES[f"{wildcards.vcf}"]['file']
    output:
        temp("tmp/variants_{vcf}/ready")
    params:
        output_prefix="tmp/variants_{vcf}/variants"
    conda: "envs/prohap.yaml"
    shell:
        script('fragment_variant_vcf.py', LEGACY_SRC) + " -i {input} -o {params.output_prefix} ; touch {output}"

rule compute_variants:
    input:
        db=REF_CHR_DB,
        tr=TRANSCRIPT_LIST,
        fasta=REF_CDNA,
        flag="tmp/variants_{vcf}/ready",
    output:
        tsv=temp("results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_chr{chr}.tsv.gz"),
        fasta=temp("results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_chr{chr}.fa" + ('.gz' if config['var_fasta_file'].endswith('.gz') else ""))
    params:
        input_vcf="tmp/variants_{vcf}/variants_chr{chr}.vcf.gz",
        output_cdna_file=("results/" + WORKING_DIR_NAME_VAR + "variants_cdna_chr{chr}.fa" + ('.gz' if config['var_cdna_file'].endswith('.gz') else "")) if (len(config['var_cdna_file']) > 0) else "",
        acc_prefix=lambda wildcards: VARIANT_VCF_FILES[f"{wildcards.vcf}"]['fasta_accession_prefix'],
        min_af=lambda wildcards: VARIANT_VCF_FILES[f"{wildcards.vcf}"]['min_af'],
        log_file="log/{vcf}_chr{chr}.log",
        #log_file="log/provar.log",
        tmp_dir="tmp/transcripts_{vcf}",
        require_start=config['var_require_start']
    conda: "envs/prohap.yaml"
    shell:
        "mkdir -p {params.tmp_dir}; mkdir -p log; mkdir -p results; " +
        script('provar.py', LEGACY_SRC) + " "
        "-i {params.input_vcf} -db {input.db} -transcripts {input.tr} -cdna {input.fasta} "
        "-chr {wildcards.chr} -acc_prefix {params.acc_prefix} -af {params.min_af} -require_start {params.require_start} "
        "-log {params.log_file} -tmp_dir {params.tmp_dir} -output_csv {output.tsv} -output_fasta {output.fasta} " + 
        ("-output_cdna_fasta \"{params.output_cdna_file}\" " if (len(config['var_cdna_file']) > 0) else "")

rule merge_var_tables_vcf:
    input:
        expand("results/" + WORKING_DIR_NAME_VAR + "/variants_{{vcf}}/variants_chr{chr}.tsv.gz", chr=CHROMOSOMES)
    output:
        temp("results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_all.tsv.gz")
    params:
        input_file_list = ','.join(expand("results/" + WORKING_DIR_NAME_VAR + "/variants_{{vcf}}/variants_chr{chr}.tsv.gz", chr=CHROMOSOMES))
    conda: "envs/prohap.yaml"
    shell:
        script('merge_tables.py', LEGACY_SRC) + " -i {params.input_file_list} -o {output}"

rule merge_var_fasta_vcf:
    input:
        expand("results/" + WORKING_DIR_NAME_VAR + "/variants_{{vcf}}/variants_chr{chr}.fa" + ('.gz' if config['var_fasta_file'].endswith('.gz') else ""), chr=CHROMOSOMES)
    output:
        temp("results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_all.fa" + ('.gz' if config['var_fasta_file'].endswith('.gz') else ""))
    params:   
        input_file_list_cdna = expand("results/" + WORKING_DIR_NAME_VAR + "/variants_cdna_chr{chr}.fa" + ('.gz' if config['var_cdna_file'].endswith('.gz') else ""), chr=CHROMOSOMES),
        output_cdna_file = "results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_cdna_all.fa" + ('.gz' if config['var_cdna_file'].endswith('.gz') else "")
    shell:
        "cat {input} > {output}" + 
        ("cat {params.input_file_list_cdna} > {params.output_cdna_file} ; rm {params.input_file_list_cdna}" if (len(config['var_cdna_file']) > 0) else "")

rule merge_var_tables:
    input:
        expand("results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_all.tsv.gz", vcf=VARIANT_VCF_FILES.keys())
    output:
        config['var_table_file']
    params:
        input_file_list = ','.join(expand("results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_all.tsv.gz", vcf=VARIANT_VCF_FILES.keys()))
    conda: "envs/prohap.yaml"
    shell:
        script('merge_tables.py', LEGACY_SRC) + " -i {params.input_file_list} -o {output}"

rule merge_var_fasta:
    input:
        expand("results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_all.fa" + ('.gz' if config['var_fasta_file'].endswith('.gz') else ""), vcf=VARIANT_VCF_FILES.keys())
    output:
        config['var_fasta_file']
    params:
        input_file_list_cdna = ','.join(expand("results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_cdna_all.fa" + ('.gz' if config['var_cdna_file'].endswith('.gz') else ""), vcf=VARIANT_VCF_FILES.keys())),
        output_cdna_file = config['var_cdna_file']
    shell:
        "cat {input} > {output}" + 
        ("cat {params.input_file_list_cdna} > {params.output_cdna_file} ; rm {params.input_file_list_cdna}" if (len(config['var_cdna_file']) > 0) else "")

rule var_fasta_remove_stop:
    input:
        config['var_fasta_file']
    output:
        temp("results/variants_all_clean.fa.gz")
    conda: "envs/prohap.yaml"
    shell:
        script('remove_stop_codons.py', LEGACY_SRC) + " -i {input} -o {output} -min_len 6 "

# ------------------------------------ ProHap rules ------------------------------------

PHASED_VCF_GLOB = (config['phased_local_path'] + config['phased_vcf_file_name']) if len(config["phased_local_path"]) > 0 else ("data/vcf/phased/" + config['phased_vcf_file_name'])

# Opt-in input preparation (F2), the same for both arms: keep only the records that overlap an exon of a listed
# transcript, read from the indexed VCF chunks. Changes the IDs given to records whose ID is '.' (line numbers).
RESTRICT_TO_EXONS = bool(config.get('restrict_to_exons', False))
if RESTRICT_TO_EXONS:
    rule restrict_to_exons:
        input:
            db=REF_CHR_DB,
            tr=TRANSCRIPT_LIST
        output:
            temp("data/vcf/phased/chr{chr}_exons.vcf.gz")
        params:
            input_file_glob=PHASED_VCF_GLOB
        threads: 2
        shell:
            "mkdir -p data/vcf/phased ; python3 " + shlex.quote(os.path.join(workflow.basedir, 'tools', 'restrict_to_exons.py')) +
            " -i \'{params.input_file_glob}\' -chr {wildcards.chr} -db {input.db} -transcripts {input.tr} -o {output} -threads {threads}"

# The chunks the filter reads; without restrict_to_exons, no rule input: the script finds them from the glob pattern.
FILTER_INPUT = ["data/vcf/phased/chr{chr}_exons.vcf.gz"] if RESTRICT_TO_EXONS else []
FILTER_GLOB = "data/vcf/phased/chr{chr}_exons.vcf.gz" if RESTRICT_TO_EXONS else PHASED_VCF_GLOB

# The legacy arm filters into a file that prohap.py reads back; the fast arm runs the filter inside compute_haplotypes
# (F3: prohap.py -filter_glob reads the filter's output through a pipe), so filter_phased_vcf runs for legacy only.
rule filter_phased_vcf:
    input:
        vcf=FILTER_INPUT
    output:
        temp("data/vcf/phased/chr{chr}_phased_filtered.vcf.gz")
    params:
        input_file_glob=FILTER_GLOB,
        AF_threshold=config['phased_min_af'],
        AF_field=config['phased_af_field']
    shell:
        "mkdir -p data/vcf/phased ; " +
        script('vcf_filter_fix.py') + " -i \'{params.input_file_glob}\' -chr {wildcards.chr} -af {params.AF_threshold} -af_field {params.AF_field} -o {output} "

rule compute_haplotypes:
    input:
        db=REF_CHR_DB,
        tr=TRANSCRIPT_LIST,
        vcf=FILTER_INPUT if IMPLEMENTATION == 'fast' else "data/vcf/phased/chr{chr}_phased_filtered.vcf.gz",
        fasta=REF_CDNA,
        samples=config['sample_metadata_file']
    output:
        csv=temp("results/" + WORKING_DIR_NAME_HAPLO + "/haplo_chr{chr}.tsv.gz"),
        fasta=temp("results/" + WORKING_DIR_NAME_HAPLO + "/haplo_chr{chr}.fa" + ('.gz' if config['haplo_fasta_file'].endswith('.gz') else "")),
    params:
        log_file="log/prohap_chr{chr}.log",
        tmp_dir="tmp/transcript_vcf_haplo",
        output_cdna_file=("results/" + WORKING_DIR_NAME_HAPLO + "/haplo_cdna_chr{chr}.fa" + ('.gz' if config['haplo_cdna_file'].endswith('.gz') else "")) if (len(config['haplo_cdna_file']) > 0) else "",
        require_start=config['haplo_require_start'],
        ignore_UTR=config['haplo_ignore_UTR'],
        skip_start_lost=config['haplo_skip_start_lost'],
        freq_threshold=config['haplo_min_freq'],
        count_threshold=config['haplo_min_count'],
        x_par1_to=config['x_par1_to'],
        x_par2_from=config['x_par2_from'],
        max_cores=config['max_cores'],
        vcf_input=("-filter_glob \'" + FILTER_GLOB + "\' -filter_af " + str(config['phased_min_af']) + " -filter_af_field "
                   + str(config['phased_af_field'])) if IMPLEMENTATION == 'fast' else "-i \"data/vcf/phased/chr{chr}_phased_filtered.vcf.gz\""
    threads: config['max_cores']
    conda: "envs/prohap.yaml"
    shell:
        "mkdir -p {params.tmp_dir}; mkdir -p log; mkdir -p results; " +
        script('prohap.py') + " "
        "{params.vcf_input} -db \"{input.db}\" -transcripts \"{input.tr}\" -cdna \"{input.fasta}\" -s \"{input.samples}\" "
        "-chr {wildcards.chr} -min_hap_freq {params.freq_threshold} -min_hap_count {params.count_threshold} "
        "-acc_prefix enshap_{wildcards.chr} -id_prefix haplo_chr{wildcards.chr} -require_start {params.require_start} -ignore_UTR {params.ignore_UTR} -skip_start_lost {params.skip_start_lost} "
        "-x_par1_to {params.x_par1_to} -x_par2_from {params.x_par2_from} -threads {params.max_cores} -log \"{params.log_file}\" -tmp_dir \"{params.tmp_dir}\" -output_csv \"{output.csv}\" -output_fasta \"{output.fasta}\" " + 
        ("-output_cdna_fasta \"{params.output_cdna_file}\" " if (len(config['haplo_cdna_file']) > 0) else "")

rule merge_haplo_tables:
    input:
        expand("results/" + WORKING_DIR_NAME_HAPLO + "/haplo_chr{chr}.tsv.gz", chr=CHROMOSOMES)
    output:
        #config['haplo_table_file']
        temp("results/" + WORKING_DIR_NAME_HAPLO + "/haplo_all.tsv.gz")
    params:
        input_file_list = ','.join(expand("results/" + WORKING_DIR_NAME_HAPLO + "/haplo_chr{chr}.tsv.gz", chr=CHROMOSOMES))
    conda: "envs/prohap.yaml"
    shell:
        script('merge_tables.py') + " -i {params.input_file_list} -o {output}"

rule extract_sample_names:
    input:
        "results/" + WORKING_DIR_NAME_HAPLO + "/haplo_all.tsv.gz"
    output:
        haplo_tsv=config['haplo_table_file'],
        samples=('.'.join(config['haplo_table_file'].split('.')[:-2]) + '_sampleIDs.tsv.gz') if config['haplo_table_file'].endswith('.gz') else ('.'.join(config['haplo_table_file'].split('.')[:-1]) + '_sampleIDs.tsv')
    conda: "envs/prohap.yaml"
    shell:
        script('haplo_extract_sample_names.py') + " -hap_tsv {input} -o {output.haplo_tsv} -samples {output.samples} "    

rule merge_fasta:
    input:
        expand(("results/" + WORKING_DIR_NAME_HAPLO + "/haplo_chr{chr}.fa" + ('.gz' if config['haplo_fasta_file'].endswith('.gz') else "")), chr=CHROMOSOMES)
    output:
        config['haplo_fasta_file']
    params:
        input_file_list_cdna = expand("results/" + WORKING_DIR_NAME_HAPLO + "/haplo_cdna_chr{chr}.fa" + ('.gz' if config['haplo_cdna_file'].endswith('.gz') else ""), chr=CHROMOSOMES),
        output_cdna_file = config['haplo_cdna_file']
    shell:
        "cat {input} > {output} ; " + 
        ("cat {params.input_file_list_cdna} > {params.output_cdna_file} ; rm {params.input_file_list_cdna}" if (len(config['haplo_cdna_file']) > 0) else "")

rule haplo_fasta_remove_stop:
    input:
        config['haplo_fasta_file']
    output:
        temp("results/haplo_all_clean.fa.gz")
    conda: "envs/prohap.yaml"
    shell:
        script('remove_stop_codons.py') + " -i {input} -o {output} -min_len 6 "

# ------------------------------------ post-processing rules ------------------------------------

rule added_fasta_remove_stop:
    input:
        config['haplo_added_fasta']
    output:
        temp("results/haplo_added_clean.fa.gz")
    conda: "envs/prohap.yaml"
    shell:
        script('remove_stop_codons.py') + " -i {input} -o {output} -min_len 6 "

rule mix_with_reference_proteome:
    input:
        in1="data/fasta/ensembl_reference_proteinDB_" + ENSEMBL_RELEASE + "_clean.fa.gz",
        in2=TAGGED_CONTAMINANTS,
        in3=expand('{proxy}', proxy=["results/variants_all_clean.fa.gz"] if config["use_ProVar"] else []),
        in4=expand('{proxy}', proxy=["results/haplo_all_clean.fa.gz"] if config["use_ProHap"] else []),
        in5=expand('{proxy}', proxy=["results/haplo_added_clean.fa.gz"] if config["add_existing_haplo"] else []),
    output:
        temp("results/ref_contam_vcf_haplo_all_clean.fa.gz")		
    run:
        shell("cat {input.in1} {input.in2} > {output}; ")
        if config["use_ProVar"]:
            shell("cat {input.in3} >> {output}")
        if config["use_ProHap"]:
            shell("cat {input.in4} >> {output}")
        if config["add_existing_haplo"]:
            shell("cat {input.in5} >> {output}")

rule merge_duplicate_seq:
    input:
        "results/ref_contam_vcf_haplo_all_clean.fa.gz"
    output:
        temp("results/ref_contam_vcf_haplo_all_nodupl.fa.gz")
        #config['final_fasta_file']                         
    conda: "envs/prohap.yaml"
    shell:
        script('merge_duplicate_seq.py') + " -i {input} -o {output} "

# UTRs in ProHap are removed by default (can be changed), but not in ProVar -> make sure all UTRs are removed
rule remove_UTR_seq:
    input:
        "results/ref_contam_vcf_haplo_all_nodupl.fa.gz"
    output:
        config['final_fasta_file']
    conda: "envs/prohap.yaml"
    shell:
        script('remove_UTR_seq.py') + " -i {input} -o {output}"

rule simpify_fasta_headers:
    input:
        fasta=config['final_fasta_file'],
        # var_table=expand('{proxy}', proxy=[config['var_table_file']] if config["use_ProVar"] else []), ENST identifiers included in the protein IDs for ProVar - no need to search for them in the TSV
        haplo_table=expand('{proxy}', proxy=[config['haplo_table_file']] if config["use_ProHap"] else ([config['haplo_added_table']] if (config['use_ProVar'] and config['add_existing_haplo']) else [])),
        annot=REF_WHOLE_DB
    output:
        fasta=('.'.join(config['final_fasta_file'].split('.')[:-2]) + '_simplified.fasta.gz') if config['final_fasta_file'].endswith('.gz') else ('.'.join(config['final_fasta_file'].split('.')[:-1]) + '_simplified.fasta'),
        header=('.'.join(config['final_fasta_file'].split('.')[:-2]) + '_header.tsv.gz') if config['final_fasta_file'].endswith('.gz') else ('.'.join(config['final_fasta_file'].split('.')[:-1]) + '_header.tsv')
    conda: "envs/prohap.yaml"
    shell:
        script('fasta_simplify_headers.py') + " -i {input.fasta} " +
        ("-hap_tsv {input.haplo_table} " if (config["use_ProHap"] or (config['use_ProVar'] and config['add_existing_haplo'])) else "") +
        "-db {input.annot} -o {output.fasta} -header {output.header}"