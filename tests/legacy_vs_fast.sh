#!/usr/bin/env bash
# Run the reference code (legacy/src), then ProHap-fast (src/), on the same data, and check that every output file
# is byte-identical. Every change to src/ must pass this.
#
#   tests/legacy_vs_fast.sh [-j CORES] [-k] [FIXTURE | TIER ...]
#
#   FIXTURE  1kg/threshold0, 1kg/template   the 1000 Genomes sample data (tests/fixtures/1kg/README.md)
#            sim/NAME                       simulated data with the config tests/fixtures/sim/NAME.yaml
#            dupseq                         merge_duplicate_seq.py and remove_UTR_seq.py on 400,000 simulated proteins
#   TIER     small    1kg/threshold0 1kg/template sim/mix-threshold0 sim/mix-template   (the default)
#            medium   sim/n2500 sim/n5000 sim/n10000 sim/whole22 dupseq
#            large    sim/n26000 sim/n100000 sim/ttn
#            all      the three tiers
#   -j       cores for each Snakemake run (default 8)
#   -k       keep the work directory even when everything is identical
#
# Measured with -j 8 on an Apple Silicon Mac: small 3 minutes, medium 28, large 65 (55 of them the 26,000-sample
# fixture). tests/fixtures/sim/README.md says what each fixture tests. For a workflow fixture, two work directories,
# legacy/ and fast/, get the same config.yaml and inputs; Snakemake runs in one with `--config
# implementation=legacy`, then in the other with `--config implementation=fast`. Then every file either run wrote
# is compared, after decompression (gzip headers carry a timestamp). Simulated data is generated on first use into
# tests/fixtures/sim/cache/ (git-ignored).
#
# Needs the ProHap environment (snakemake on PATH) and the pre-built references of Ensembl release 113 in
# $PROHAP_REFS/ensembl.v113 (e.g. PROHAP_REFS=~/prohap-refs); no container. Exit status 0 when every output is
# identical, 1 when one differs or a run fails; the work directory is then kept. A fast run in which a script handed
# its input over to the reference code (src/modules/fallback.py) counts as failed: it would compare the reference
# with itself.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
FALLBACK='ProHap-fast: falling back to the reference code'     # src/modules/fallback.py
SMALL="1kg/threshold0 1kg/template sim/mix-threshold0 sim/mix-template"
MEDIUM="sim/n2500 sim/n5000 sim/n10000 sim/whole22 dupseq"
LARGE="sim/n26000 sim/n100000 sim/ttn"

cores=8
keep=0
while getopts j:kh opt; do
    case $opt in
        j) cores=$OPTARG ;;
        k) keep=1 ;;
        *) sed -n '2,28p' "$0"; exit 2 ;;
    esac
done
shift $((OPTIND - 1))
[ $# -gt 0 ] || set -- small
fixtures=()
for arg in "$@"; do
    case $arg in
        small) fixtures+=($SMALL) ;;
        medium) fixtures+=($MEDIUM) ;;
        large) fixtures+=($LARGE) ;;
        all) fixtures+=($SMALL $MEDIUM $LARGE) ;;
        *) fixtures+=("$arg") ;;
    esac
done
for fixture in "${fixtures[@]}"; do
    case $fixture in
        1kg/* | sim/*) [ -f "$ROOT/tests/fixtures/$fixture.yaml" ] && continue ;;
        dupseq) continue ;;
    esac
    echo "no fixture $fixture (tests/fixtures/$fixture.yaml)" >&2
    exit 2
done
command -v snakemake > /dev/null || { echo 'snakemake is not on PATH: activate the ProHap environment' >&2; exit 2; }

work=$(mktemp -d "${TMPDIR:-/tmp}/legacy_vs_fast.XXXXXX")
echo "Work directory: $work"

# The 1kG inputs, shared by its runs: the VCFs of chromosomes 1 (two chunks), 21, 22 and X, the sample metadata, the
# contaminants and the transcript list.
input_1kg() {
    if [ ! -d "$work/input/1kg" ]; then
        mkdir -p "$work/input/1kg/sample_data"
        tar -xzf "$ROOT/sample_data/sample_1kGP_common_global.tar.gz" -C "$work/input/1kg/sample_data" \
            sample_chr1_1kGP_GRCh38_phased.vcf sample_chr1_1kGP_GRCh38_phased_2.vcf \
            sample_chr21_1kGP_GRCh38_phased.vcf sample_chr22_1kGP_GRCh38_phased.vcf \
            sample_chrX_1kGP_GRCh38_phased.vcf
        cp "$ROOT/igsr_samples.tsv" "$ROOT/crap.fasta" "$ROOT/tests/fixtures/1kg/transcripts.csv" "$work/input/1kg/"
    fi
}

# The directory of a simulated dataset, generated into the cache unless it is there already.
sim_data() {
    python3 "$ROOT/tests/fixtures/sim/make_data.py" "$1"
}

# Every file a run wrote: everything in its directory except its config, the links to the inputs and Snakemake's
# own records in .snakemake/.
outputs() {
    (cd "$1" && find . -type f ! -path ./config.yaml ! -path './.snakemake/*' | sed 's|^\./||')
}

is_gzip() {
    [ "$(head -c 2 "$1" | od -An -tx1 | tr -d ' \n')" = 1f8b ]
}

# Show the end of a failed run's log: the job's own error comes just before Snakemake's "Error in rule"; a config
# error is at the end.
show_failure() {
    local line
    line=$(grep -n -m 1 '^Error in rule' "$1" | cut -d: -f1 || true)
    if [ -n "$line" ]; then
        sed -n "$((line > 15 ? line - 15 : 1)),$((line + 1))p" "$1"
    else
        tail -n 20 "$1"
    fi | cut -c 1-300 | sed 's/^/    /'
}

# run_arm FIXTURE ARM DIR: one run of FIXTURE with the code of ARM, in DIR (which holds its inputs); its output goes
# to DIR.log.
run_arm() {
    local fixture=$1 arm=$2 dir=$3 src
    if [ "$fixture" = dupseq ]; then
        src=$ROOT/src
        [ "$arm" = fast ] || src=$ROOT/legacy/src
        (cd "$dir" && export PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 &&
            python3 "$src/merge_duplicate_seq.py" -i input.fa.gz -o nodupl.fa.gz &&
            python3 "$src/remove_UTR_seq.py" -i nodupl.fa.gz -o final.fa)
    else
        PYTHONHASHSEED=0 snakemake --snakefile "$ROOT/Snakefile" --directory "$dir" --cores "$cores" --notemp \
            --config implementation="$arm"
    fi
}

failed=0
total=0
different=0
start_all=$SECONDS
for fixture in "${fixtures[@]}"; do
    echo
    echo "== $fixture"
    case $fixture in
        1kg/*) input_1kg ;;
        # the dataset of sim/NAME-VARIANT.yaml is NAME
        sim/*) name=${fixture#sim/}; data=$(sim_data "${name%%-*}") || exit 2 ;;
        dupseq) data=$(sim_data dupseq) || exit 2 ;;
    esac
    fdir=$work/${fixture//\//_}
    for arm in legacy fast; do
        dir=$fdir/$arm
        mkdir -p "$dir"
        case $fixture in
            1kg/*)
                cp "$ROOT/tests/fixtures/$fixture.yaml" "$dir/config.yaml"
                ln -s "$work"/input/1kg/{sample_data,igsr_samples.tsv,crap.fasta,transcripts.csv} "$dir/" ;;
            sim/*)
                cp "$ROOT/tests/fixtures/$fixture.yaml" "$dir/config.yaml"
                ln -s "$data" "$dir/simulated"
                ln -s "$ROOT/crap.fasta" "$dir/crap.fasta" ;;
            dupseq)
                ln -s "$data/input.fa.gz" "$dir/input.fa.gz" ;;
        esac
        start=$SECONDS
        if ! run_arm "$fixture" "$arm" "$dir" > "$fdir/$arm.log" 2>&1; then
            echo "$arm run FAILED (log: $fdir/$arm.log):"
            show_failure "$fdir/$arm.log"
            failed=1
            continue 2
        fi
        printf '%-7s run: %4d s\n' "$arm" $((SECONDS - start))
        if [ "$arm" = fast ] && grep -q "$FALLBACK" "$fdir/$arm.log"; then
            echo "fast run FELL BACK to the reference code (log: $fdir/$arm.log):"
            grep "$FALLBACK" "$fdir/$arm.log" | cut -c 1-300 | sed 's/^/    /'
            failed=1
        fi
    done

    legacy=$fdir/legacy
    fast=$fdir/fast
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
elapsed=$((SECONDS - start_all))
if [ $failed = 0 ] && [ $different = 0 ]; then
    printf 'IDENTICAL: all %d outputs of %d fixture(s), in %d min %d s.\n' $total ${#fixtures[@]} \
        $((elapsed / 60)) $((elapsed % 60))
    if [ $keep = 0 ]; then rm -rf "$work"; else echo "Kept $work"; fi
    exit 0
fi
[ $failed = 0 ] || echo 'A run FAILED or fell back to the reference code.'
[ $different = 0 ] || echo "DIFFERENT: $different of $total outputs."
echo "Kept $work"
exit 1
