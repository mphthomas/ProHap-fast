#!/usr/bin/env bash
#
# ProHap container entrypoint.
#
# Runs the workflow from the read-only copy of the repository baked into the
# image (/opt/prohap) against a writable working directory supplied by the
# caller. Nothing belonging to the image is ever linked or copied into that
# directory, so every file the run produces is a real file readable on the host
# after the container exits.
#
# The working directory is simply the current directory, which under
# Singularity is your host directory, already bind-mounted. So:
#
#   cd /scratch/prohap_run        # contains config.yaml
#   singularity run --bind /data/phased_vcfs IMAGE snakemake --cores 16 -p
#
# Input VCFs and outputs live wherever config.yaml says; bind those paths in
# and refer to them by their normal host paths.
#
# Environment (all optional):
#   PROHAP_WORK      working directory            (default: current directory)
#   PROHAP_REFS      gzipped Ensembl references   (default /opt/prohap/refs)
#   PROHAP_HOME      repository inside the image  (default /opt/prohap)
#
# Any arguments are passed through to snakemake. With none, a dry run is
# printed rather than silently starting a multi-day job.

set -euo pipefail

PROHAP_HOME="${PROHAP_HOME:-/opt/prohap}"
PROHAP_WORK="${PROHAP_WORK:-$PWD}"
PROHAP_REFS="${PROHAP_REFS:-${PROHAP_HOME}/refs}"

# Singularity does not reliably propagate the image's ENV, so do not depend on
# PATH having survived the transition from Docker.
export PATH="/opt/conda/envs/ProHap/bin:${PATH:-/usr/local/bin:/usr/bin:/bin}"
# pandas and pyarrow are pip wheels: their C++ extensions must load the environment's libstdc++, not the image's
# older system copy, or on linux/arm64 importing gffutils (sqlite3, then ICU) after pandas fails.
export LD_LIBRARY_PATH="/opt/conda/envs/ProHap/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

die() { echo "[prohap] ERROR: $*" >&2; exit 1; }
log() { echo "[prohap] $*"; }

# Escape hatch for debugging the image itself.
if [ "${1:-}" = "shell" ]; then
    shift
    exec bash "$@"
fi

# `... IMAGE snakemake --cores 16 -p` is the natural way to write this and was
# what the previous entrypoint required, so accept it and drop the redundant
# leading token. The remaining arguments are snakemake's.
if [ "${1:-}" = "snakemake" ]; then
    shift
fi

command -v snakemake >/dev/null 2>&1 \
    || die "snakemake not found on PATH (${PATH})"

[ -d "$PROHAP_WORK" ] \
    || die "work directory '$PROHAP_WORK' does not exist. Bind-mount it, or set PROHAP_WORK."
[ -w "$PROHAP_WORK" ] \
    || die "work directory '$PROHAP_WORK' is not writable by uid $(id -u). Check the mount and ownership."

# The Snakefile carries `configfile: "config.yaml"`, resolved against the
# working directory, so that is where the config lives. Deliberately not
# copied in from elsewhere: writing over a file called config.yaml in the
# caller's own directory is too surprising to be worth the flexibility.
PROHAP_CONFIG="${PROHAP_WORK}/config.yaml"
[ -f "$PROHAP_CONFIG" ] \
    || die "no config.yaml in the work directory '$PROHAP_WORK'.
       Copy config_tre_template.yaml there as config.yaml and edit it."
log "using config $PROHAP_CONFIG"

# Snakemake and conda want somewhere to write; under `--user` on the host, or
# under Singularity with --no-home, the real HOME may not be writable.
if [ ! -w "${HOME:-/nonexistent}" ]; then
    export HOME="${PROHAP_WORK}/.home"
    mkdir -p "$HOME"
fi
export XDG_CACHE_HOME="${PROHAP_WORK}/.cache"
mkdir -p "$XDG_CACHE_HOME"

# Put the reference files where the download_* rules would have written them,
# so those rules never fire and never need the network.
python3 "${PROHAP_HOME}/container/stage_refs.py" \
    --config "$PROHAP_CONFIG" \
    --refs   "$PROHAP_REFS" \
    --workdir "$PROHAP_WORK"

if [ "$#" -eq 0 ]; then
    log "no snakemake arguments given - showing a dry run."
    log "re-run with e.g. '--cores 16 -p' to execute."
    set -- --dry-run
fi

log "work directory: $PROHAP_WORK"
log "snakemake $*"

# --rerun-triggers mtime: without this, editing the workflow invalidates the
# cached gffutils databases and forces a multi-hour reparse of the GTF on an
# image upgrade. Listed before "$@" so an explicit override still wins.
exec snakemake \
    --snakefile "${PROHAP_HOME}/Snakefile" \
    --directory "$PROHAP_WORK" \
    --rerun-triggers mtime \
    "$@"
