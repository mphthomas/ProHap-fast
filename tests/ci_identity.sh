#!/usr/bin/env bash
# The identity checks of CI (.github/workflows/identity.yml), run inside a prohap-fast image under Docker or
# Apptainer. Also runs on a workstation with Docker (Docker Desktop on a Mac runs the linux/amd64 image with Rosetta).
#
#   tests/ci_identity.sh --runtime docker|apptainer --image IMAGE --work DIR verify
#   tests/ci_identity.sh --runtime docker|apptainer --image IMAGE --work DIR tiers TIER... [--write-golden]
#
#   verify   `verify` through the image's entrypoint, as a user runs it (this user, the work directory mounted at its
#            own path): the 1kG fixture tests/fixtures/1kg/threshold0.yaml on chromosomes 1 (two chunks) and 22,
#            --samples 300 --transcripts 30. Every output must be IDENTICAL.
#   tiers    tests/legacy_vs_fast.sh -g TIER... with the image's code, environment and references
#            (PROHAP_CODE=/opt/prohap; the fixtures and this script come from this checkout, mounted): every output
#            IDENTICAL, and the reference's outputs equal to the golden checksums tests/fixtures/expected.v113.sha256.
#            --write-golden writes them instead (-G): only for a deliberate change, an output change with its
#            CHANGELOG entry, a new Ensembl release or a new fixture, and only from the linux/amd64 image.
#
# IMAGE is a Docker image reference, or a .sif file for Apptainer. DIR holds the work directories (verify/, tmp/ with
# legacy_vs_fast.sh's work directory) and the logs. Exit status 0 when everything is identical.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
runtime='' image='' work='' write_golden=0
while [ $# -gt 0 ]; do
    case $1 in
        --runtime) runtime=$2; shift 2 ;;
        --image) image=$2; shift 2 ;;
        --work) work=$2; shift 2 ;;
        --write-golden) write_golden=1; shift ;;
        -h | --help) sed -n '2,/^set -e/p' "$0" | sed '$d; s/^# \{0,1\}//'; exit 0 ;;
        *) break ;;
    esac
done
[ $# -gt 0 ] || { echo 'ci_identity.sh: verify or tiers TIER... (see --help)' >&2; exit 2; }
case $runtime in docker | apptainer) ;; *) echo 'ci_identity.sh: --runtime docker or apptainer' >&2; exit 2 ;; esac
[ -n "$image" ] && [ -n "$work" ] || { echo 'ci_identity.sh: --image and --work are required' >&2; exit 2; }
mkdir -p "$work"
work=$(cd "$work" && pwd)

# run_in_image DIR ARGS...: the image's entrypoint with ARGS, in DIR, as this user; DIR, this checkout and the work
# directory mounted at their own paths
run_in_image() {
    local dir=$1
    shift
    if [ "$runtime" = docker ]; then
        docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -e TMPDIR="$work/tmp" -e PROHAP_CODE=/opt/prohap \
            -v "$ROOT:$ROOT" -v "$work:$work" -w "$dir" "$image" "$@"
    else
        apptainer run --bind "$ROOT" --bind "$work" --pwd "$dir" \
            --env "TMPDIR=$work/tmp,PROHAP_CODE=/opt/prohap" "$image" "$@"
    fi
}

case $1 in
    verify)
        dir=$work/verify
        rm -rf "$dir"
        mkdir -p "$dir/sample_data"
        tar -xzf "$ROOT/sample_data/sample_1kGP_common_global.tar.gz" -C "$dir/sample_data" \
            sample_chr1_1kGP_GRCh38_phased.vcf sample_chr1_1kGP_GRCh38_phased_2.vcf \
            sample_chr22_1kGP_GRCh38_phased.vcf
        cp "$ROOT/igsr_samples.tsv" "$ROOT/crap.fasta" "$ROOT/tests/fixtures/1kg/transcripts.csv" "$dir/"
        cp "$ROOT/tests/fixtures/1kg/threshold0.yaml" "$dir/config.yaml"
        run_in_image "$dir" verify --cores 4 --chromosomes 1,22 --samples 300 --transcripts 30 2>&1 \
            | tee "$work/verify.log"
        ;;
    tiers)
        shift
        mkdir -p "$work/tmp"
        flag=-g
        [ $write_golden = 0 ] || flag=-G
        run_in_image "$ROOT" shell "$ROOT/tests/legacy_vs_fast.sh" -j 4 $flag "$@" 2>&1 | tee "$work/tiers.log"
        ;;
    *) echo "ci_identity.sh: unknown command $1 (verify or tiers)" >&2; exit 2 ;;
esac
