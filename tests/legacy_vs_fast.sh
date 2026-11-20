#!/usr/bin/env bash
# Run the reference code (legacy/src), then ProHap-fast (src/), on the same data, and check that every output file
# is byte-identical. Every change to src/ must pass this.
#
#   tests/legacy_vs_fast.sh [-j CORES] [-k] [-g | -G] [FIXTURE | TIER ...]
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
#   -g       also check the reference's outputs against the golden checksums, tests/fixtures/expected.v113.sha256
#            (Linux only). Fast vs legacy cannot see a change that hits both arms (the Snakefile, the references,
#            the environment); the golden checksums can. CI runs this inside the image (tests/ci_identity.sh).
#   -G       write the golden checksums of the fixtures that ran instead, keeping the others' (Linux only; only for
#            a deliberate change: an output change with its CHANGELOG entry, a new release or a new fixture)
#
# Measured with -j 8 on an Apple Silicon Mac (0.2.0): small 4 minutes, medium 15, large 34 (28 of them the reference
# on the 26,000-sample fixture); a fast run takes 15-35 s, mostly steps both arms share. tests/fixtures/sim/README.md
# says what each fixture tests. For a workflow fixture, two work directories,
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
#
# PROHAP_CODE, if set, is the directory whose Snakefile, src/ and legacy/src run (CI sets /opt/prohap, the image's
# copy); the fixtures always come from this checkout.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
CODE=${PROHAP_CODE:-$ROOT}
GOLDEN=$ROOT/tests/fixtures/expected.v113.sha256
FALLBACK='ProHap-fast: falling back to the reference code'     # src/modules/fallback.py
SMALL="1kg/threshold0 1kg/template sim/mix-threshold0 sim/mix-template"
MEDIUM="sim/n2500 sim/n5000 sim/n10000 sim/whole22 dupseq"
LARGE="sim/n26000 sim/n100000 sim/ttn"

cores=8
keep=0
golden=''
while getopts j:kgGh opt; do
    case $opt in
        j) cores=$OPTARG ;;
        k) keep=1 ;;
        g) golden=check ;;
        G) golden=write ;;
        *) sed -n '2,/^set -euo/p' "$0" | sed '$d'; exit 2 ;;
    esac
done
if [ -n "$golden" ] && [ "$(uname -s)" != Linux ]; then
    echo "-g and -G need Linux: the golden checksums come from the Linux images, and the macOS build of pandas" \
         "parses some floats differently" >&2
    exit 2
fi
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

# Every file a run wrote: everything in its directory except its config, the links to the inputs, Snakemake's
# own records in .snakemake/, and the filtered VCF, which only the legacy arm writes (the fast arm filters inside
# prohap.py, F3).
outputs() {
    (cd "$1" && find . -type f ! -path ./config.yaml ! -path './.snakemake/*' \
        ! -path './data/vcf/phased/*_phased_filtered.vcf.gz' | sed 's|^\./||')
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
        src=$CODE/src
        [ "$arm" = fast ] || src=$CODE/legacy/src
        (cd "$dir" && export PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 &&
            python3 "$src/merge_duplicate_seq.py" -i input.fa.gz -o nodupl.fa.gz &&
            python3 "$src/remove_UTR_seq.py" -i nodupl.fa.gz -o final.fa)
    else
        PYTHONHASHSEED=0 snakemake --snakefile "$CODE/Snakefile" --directory "$dir" --cores "$cores" --notemp \
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
        if [ -n "$golden" ] && [ -f "$legacy/$file" ]; then
            printf '%s  %s\n' "$(gzip -cdf "$legacy/$file" | sha256sum | cut -d ' ' -f 1)" "$fixture/$file" \
                >> "$work/legacy.sha256"
        fi
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

# The golden checksums: tests/fixtures/expected.v113.sha256 holds "<sha256>  <fixture>/<file>" for the reference's
# outputs; only the entries of the fixtures that ran are checked or written.
golden_failed=0
if [ -n "$golden" ] && [ $failed = 0 ] && [ $different = 0 ]; then
    echo
    if ! python3 - "$golden" "$GOLDEN" "$work/legacy.sha256" "${fixtures[@]}" <<'PYTHON'
import os
import sys

mode, path, run_path, fixtures = sys.argv[1], sys.argv[2], sys.argv[3], set(sys.argv[4:])


def read(p):
    sums = {}
    if os.path.exists(p):
        for line in open(p):
            if line.strip() and not line.startswith('#'):
                digest, name = line.rstrip('\n').split('  ', 1)
                sums[name] = digest
    return sums


ran = read(run_path)
golden = read(path)
mine = lambda name: any(name.startswith(f + '/') for f in fixtures)
if mode == 'write':
    merged = {k: v for k, v in golden.items() if not mine(k)}
    merged.update(ran)
    with open(path, 'w') as fh:
        fh.write('# sha256 of the decompressed outputs of the reference (legacy arm) on the fixtures of\n'
                 '# tests/legacy_vs_fast.sh, Ensembl release 113, written in the linux/amd64 image by\n'
                 '# tests/ci_identity.sh --write-golden (legacy_vs_fast.sh -G). Checked by CI (-g). Linux only.\n')
        for name in sorted(merged):
            fh.write('%s  %s\n' % (merged[name], name))
    print('Wrote %d golden checksums (%d of the fixtures that ran) to %s' % (len(merged), len(ran), path))
    sys.exit(0)
expected = {k: v for k, v in golden.items() if mine(k)}
problems = ['%s: in the golden checksums, not written by this run' % k for k in sorted(set(expected) - set(ran))]
problems += ['%s: written by this run, not in the golden checksums' % k for k in sorted(set(ran) - set(expected))]
problems += ['%s: sha256 %s, golden %s' % (k, ran[k][:12], expected[k][:12])
             for k in sorted(set(ran) & set(expected)) if ran[k] != expected[k]]
if problems:
    print('GOLDEN: the reference\'s outputs differ from %s:' % path)
    for line in problems:
        print('  ' + line)
    sys.exit(1)
print('GOLDEN: all %d outputs of the reference match the golden checksums.' % len(ran))
PYTHON
    then
        golden_failed=1
    fi
fi

echo
elapsed=$((SECONDS - start_all))
if [ $failed = 0 ] && [ $different = 0 ] && [ $golden_failed = 0 ]; then
    printf 'IDENTICAL: all %d outputs of %d fixture(s), in %d min %d s.\n' $total ${#fixtures[@]} \
        $((elapsed / 60)) $((elapsed % 60))
    if [ $keep = 0 ]; then rm -rf "$work"; else echo "Kept $work"; fi
    exit 0
fi
[ $failed = 0 ] || echo 'A run FAILED or fell back to the reference code.'
[ $different = 0 ] || echo "DIFFERENT: $different of $total outputs."
[ $golden_failed = 0 ] || echo "Every output IDENTICAL, but the reference's outputs differ from the golden checksums."
echo "Kept $work"
exit 1
