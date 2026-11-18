#!/usr/bin/env bash
# Run the reference code (legacy/src), then ProHap-fast (src/), on the same data, and check that every output file
# is byte-identical. Every change to src/ must pass this.
#
#   tests/legacy_vs_fast.sh [-j CORES] [-k] [CONFIG ...]
#
#   CONFIG  fixture configs in tests/fixtures/1kg/ to run (default: threshold0 template)
#   -j      cores for each Snakemake run (default 8)
#   -k      keep the work directory even when everything is identical
#
# The data is the 1kG fixture (tests/fixtures/1kg/README.md). For each config, two work directories, legacy/ and
# fast/, get the same config.yaml and inputs; Snakemake runs in one with `--config implementation=legacy`, then in
# the other with `--config implementation=fast`. Then every file either run wrote is compared, after decompression
# (gzip headers carry a timestamp). About a minute and a half on 8 cores.
#
# Needs the ProHap environment (snakemake on PATH) and the pre-built references of Ensembl release 113 in
# $PROHAP_REFS/ensembl.v113 (e.g. PROHAP_REFS=~/prohap-refs); no container. Exit status 0 when every output is
# identical, 1 when one differs or a run fails; the work directory is then kept.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
FIXTURE=$ROOT/tests/fixtures/1kg

cores=8
keep=0
while getopts j:kh opt; do
    case $opt in
        j) cores=$OPTARG ;;
        k) keep=1 ;;
        *) sed -n '2,18p' "$0"; exit 2 ;;
    esac
done
shift $((OPTIND - 1))
[ $# -gt 0 ] || set -- threshold0 template
for config in "$@"; do
    [ -f "$FIXTURE/$config.yaml" ] || { echo "no config $FIXTURE/$config.yaml" >&2; exit 2; }
done
command -v snakemake > /dev/null || { echo 'snakemake is not on PATH: activate the ProHap environment' >&2; exit 2; }

work=$(mktemp -d "${TMPDIR:-/tmp}/legacy_vs_fast.XXXXXX")
echo "Work directory: $work"

# The inputs, shared by every run: the VCFs of chromosomes 1 (two chunks), 21, 22 and X, the sample metadata, the
# contaminants and the transcript list.
mkdir -p "$work/input/sample_data"
tar -xzf "$ROOT/sample_data/sample_1kGP_common_global.tar.gz" -C "$work/input/sample_data" \
    sample_chr1_1kGP_GRCh38_phased.vcf sample_chr1_1kGP_GRCh38_phased_2.vcf sample_chr21_1kGP_GRCh38_phased.vcf \
    sample_chr22_1kGP_GRCh38_phased.vcf sample_chrX_1kGP_GRCh38_phased.vcf
cp "$ROOT/igsr_samples.tsv" "$ROOT/crap.fasta" "$FIXTURE/transcripts.csv" "$work/input/"

# Every file a run wrote: everything in its work directory except its config, the input links and Snakemake's own
# records in .snakemake/.
outputs() {
    (cd "$1" && find . -type f ! -path ./config.yaml ! -path './.snakemake/*' | sed 's|^\./||')
}

is_gzip() {
    [ "$(head -c 2 "$1" | od -An -tx1 | tr -d ' \n')" = 1f8b ]
}

failed=0
total=0
different=0
for config in "$@"; do
    echo
    echo "== $config"
    for arm in legacy fast; do
        dir=$work/$config/$arm
        mkdir -p "$dir"
        cp "$FIXTURE/$config.yaml" "$dir/config.yaml"
        ln -s ../../input/sample_data ../../input/igsr_samples.tsv ../../input/crap.fasta \
            ../../input/transcripts.csv "$dir/"
        log=$work/$config/$arm.log
        start=$SECONDS
        if ! PYTHONHASHSEED=0 snakemake --snakefile "$ROOT/Snakefile" --directory "$dir" --cores "$cores" \
                --notemp --config implementation=$arm > "$log" 2>&1; then
            # the job's own error comes just before Snakemake's "Error in rule"; a config error is at the end
            echo "$arm run FAILED (log: $log):"
            line=$(grep -n -m 1 '^Error in rule' "$log" | cut -d: -f1 || true)
            if [ -n "$line" ]; then
                sed -n "$((line > 15 ? line - 15 : 1)),$((line + 1))p" "$log"
            else
                tail -n 20 "$log"
            fi | cut -c 1-300 | sed 's/^/    /'
            failed=1
            continue 2
        fi
        printf '%-7s run: %4d s\n' "$arm" $((SECONDS - start))
    done

    legacy=$work/$config/legacy
    fast=$work/$config/fast
    for file in $( (outputs "$legacy"; outputs "$fast") | sort -u); do
        total=$((total + 1))
        if [ ! -f "$fast/$file" ]; then
            result='only in legacy'
        elif [ ! -f "$legacy/$file" ]; then
            result='only in fast'
        elif is_gzip "$legacy/$file" && ! is_gzip "$fast/$file"; then
            result='DIFFERENT (gzip in legacy only)'
        elif ! is_gzip "$legacy/$file" && is_gzip "$fast/$file"; then
            result='DIFFERENT (gzip in fast only)'
        elif cmp -s <(gzip -cdf "$legacy/$file") <(gzip -cdf "$fast/$file"); then
            result=IDENTICAL
        else
            result=DIFFERENT
        fi
        printf '  %-10s %s\n' "$result" "$file"
        if [ "$result" != IDENTICAL ]; then
            different=$((different + 1))
            if [ "$result" = DIFFERENT ]; then
                # the first difference (< legacy, > fast), long lines cut at 160 characters
                diff <(gzip -cdf "$legacy/$file") <(gzip -cdf "$fast/$file") | head -n 4 | cut -c 1-160 \
                    | sed 's/^/        /' || true
            fi
        fi
    done
done

echo
if [ $failed = 0 ] && [ $different = 0 ]; then
    echo "IDENTICAL: all $total outputs of $# config(s)."
    if [ $keep = 0 ]; then rm -rf "$work"; else echo "Kept $work"; fi
    exit 0
fi
[ $failed = 0 ] || echo 'A run FAILED.'
[ $different = 0 ] || echo "DIFFERENT: $different of $total outputs."
echo "Kept $work"
exit 1
