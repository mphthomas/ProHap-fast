#!/usr/bin/env bash
#
# Simple script to compare legacy and fast implementations of ProHap on simulated datasets.
# Run both, compare outputs are byte identical, remove the temporary files.
#
# This can also be done directly with snakemake:
#   snakemake --snakefile $SNAKEFILE --configfile $CONFIG --directory run_legacy --cores 4 -p --config implementation=legacy
#   snakemake --snakefile $SNAKEFILE --configfile $CONFIG --directory run_fast   --cores 4 -p --config implementation=fast

# But to avoid any ambiguity this script directly runs both implementations using the same arguments.
# To avoid downloading dependencies, run this with `docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":"$PWD" -w "$PWD" $IMAGE shell ./legacy_vs_fast.sh`.
# On my local machine I set IMAGE='mphthomas/prohap-fast:1.0.0-ensembl.v113-arm64'


export PYTHONHASHSEED=0
# Runs for ensembl version 113 on chromosome 22
REF=${PROHAP_REFS:-$HOME/prohap-refs}/ensembl.v113
OUT="legacy_vs_fast_out"
CHR=22
METADATA=igsr_samples.tsv
INPUT_DATA='sample_data/sample_chr22_1kGP_GRCh38_phased*.vcf'
mkdir -p legacy_vs_fast_out/out_legacy/tmp legacy_vs_fast_out/out_fast/tmp

# These are usually specified in config but defined here for repeated use in this script.
prohap_args() {   # $1 = output directory given by prohap_args below
  ARGS=(-i "$1/filtered_chr${CHR}.vcf.gz"
        -db "$REF/Homo_sapiens.GRCh38.113.chr_patch_hapl_scaff_chr${CHR}.db"
        -transcripts "$REF/transcripts_reference_113.csv" -cdna "$REF/total_cdnas_113.fa" -s $METADATA
        -chr ${CHR} -min_hap_freq 0 -min_hap_count 0 -acc_prefix enshap_${CHR} -id_prefix haplo_chr${CHR}
        -require_start 1 -ignore_UTR 1 -skip_start_lost 1 -x_par1_to 2781479 -x_par2_from 155701383
        -threads 1 -log "$1/prohap_chr${CHR}.log" -tmp_dir "$1/tmp"
        -output_csv "$1/haplo_chr${CHR}.tsv.gz" -output_fasta "$1/haplo_chr${CHR}.fa")
}

### Comparison, only change is the legacy in path name.
{ 
    # Legacy
    start=`date +%s`
    python3 tools/run_fork.py legacy/src/vcf_filter_fix.py \
    -i "$INPUT_DATA" -chr $CHR -af 0.01 -af_field AF \
    -o "$OUT/out_legacy/filtered_chr${CHR}.vcf.gz"
    prohap_args "$OUT/out_legacy"; python3 tools/run_fork.py legacy/src/prohap.py "${ARGS[@]}"
    end=`date +%s`
    legacy_runtime=$((end-start))

    # Fast
    start=`date +%s`
    python3 tools/run_fork.py src/vcf_filter_fix.py \
    -i "$INPUT_DATA" -chr $CHR -af 0.01 -af_field AF \
    -o "$OUT/out_fast/filtered_chr${CHR}.vcf.gz"
    prohap_args "$OUT/out_fast";   python3 tools/run_fork.py src/prohap.py "${ARGS[@]}"
    end=`date +%s`
    fast_runtime=$((end-start))
    
    # Compare outputs
    for f in filtered_chr${CHR}.vcf.gz haplo_chr${CHR}.tsv.gz haplo_chr${CHR}.fa prohap_chr${CHR}.log; do
    cmp <(gzip -dcf "$OUT/out_legacy/$f") <(gzip -dcf "$OUT/out_fast/$f") && echo "$f identical"
    done

    echo "Legacy runtime: $legacy_runtime second(s)"
    echo "Fast runtime: $fast_runtime second(s)";
    exit; 
}

