#!/usr/bin/env bash
#
# prohap-fast container entrypoint.
#
# The work directory is the current directory. It holds config.yaml, the intermediates and the logs; outputs go
# wherever config.yaml says. The workflow, its environment and the pre-built references stay in the image and are
# read in place: nothing is linked or copied into the work directory, and nothing is downloaded.
#
#   cd /path/to/workdir                                   # contains config.yaml
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":"$PWD" -w "$PWD" -v /data:/data:ro IMAGE [COMMAND] [ARGS]
#   singularity run --bind /data IMAGE [COMMAND] [ARGS]  # the current directory and user come automatically
#
# Mount every path config.yaml refers to at the same path inside the container.
#
# COMMAND:
#   (none)              a dry run: what the workflow would do, without doing it
#   [snakemake] ARGS    run the workflow: snakemake ARGS, e.g. --cores 16 -p
#   verify [ARGS]       run the reference and ProHap-fast on the same input and compare every output
#                       (/opt/prohap/verify/README.md; verify --help)
#   shell [ARGS]        bash, for looking around inside the image
#   help                this text
#
# PROHAP_WORK, if set, is the work directory instead of the current directory.

set -euo pipefail

PROHAP_HOME="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
PROHAP_WORK="${PROHAP_WORK:-$PWD}"

# Singularity passes the host's environment in and does not reliably apply the image's ENV, so everything a run
# depends on is set here. The references are always the image's (a config's refs_dir can still point elsewhere):
# a PROHAP_REFS from the host must not leak in. PYTHONNOUSERSITE keeps a ~/.local site-packages, which Singularity
# mounts with $HOME, from shadowing the environment's packages.
export PATH="/opt/conda/envs/ProHap/bin:${PATH:-/usr/local/bin:/usr/bin:/bin}"
# pandas and pyarrow are pip wheels: their C++ extensions must load the environment's libstdc++, not the image's
# older system copy, or on linux/arm64 importing gffutils (sqlite3, then ICU) after pandas fails.
export LD_LIBRARY_PATH="/opt/conda/envs/ProHap/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export PROHAP_HOME
export PROHAP_REFS="${PROHAP_HOME}/refs"
export PYTHONHASHSEED=0
export PYTHONNOUSERSITE=1

die() { echo "[prohap] ERROR: $*" >&2; exit 1; }
log() { echo "[prohap] $*" >&2; }
usage() { sed -n '3,/^set -euo/p' "${BASH_SOURCE[0]}" | sed '$d; s/^# \{0,1\}//'; }

command="${1:-}"
case "$command" in
    shell) shift; exec bash "$@" ;;
    help | -h | --help) usage; exit 0 ;;
    verify) shift ;;
    snakemake) shift; command=run ;;
    *) command=run ;;
esac

# The release of the references in this image: the one ensembl.v<release>/ directory under refs/.
releases=("$PROHAP_REFS"/ensembl.v*/release.txt)
[ "${#releases[@]}" -eq 1 ] && [ -f "${releases[0]}" ] \
    || die "the image should hold the references of exactly one Ensembl release in $PROHAP_REFS"
export ENSEMBL_RELEASE="$(cat "${releases[0]}")"
version="$(cat "$PROHAP_HOME/VERSION")"
image="prohap-fast:${version}-ensembl.v${ENSEMBL_RELEASE}"

[ -d "$PROHAP_WORK" ] \
    || die "the work directory $PROHAP_WORK does not exist. Mount it (docker run -v \"\$PWD\":\"\$PWD\" -w \"\$PWD\"), or set PROHAP_WORK."
[ -w "$PROHAP_WORK" ] \
    || die "the work directory $PROHAP_WORK is not writable by uid $(id -u). Check the mount, and run Docker with -u \"\$(id -u):\$(id -g)\"."
cd "$PROHAP_WORK"

# The Snakefile and verify read config.yaml from the work directory. It is never written or copied from elsewhere.
[ -f config.yaml ] \
    || die "no config.yaml in the work directory $PROHAP_WORK.
       Start from the template in the image: $PROHAP_HOME/config_tre_template.yaml, e.g.
       docker run --rm -u \"\$(id -u):\$(id -g)\" -v \"\$PWD\":\"\$PWD\" -w \"\$PWD\" $image shell -c 'cp $PROHAP_HOME/config_tre_template.yaml config.yaml'"
python3 "$PROHAP_HOME/container/check_release.py" config.yaml "$ENSEMBL_RELEASE" "$image"

# Snakemake wants somewhere to write. Under docker run -u, or Singularity with --no-home, HOME may not be writable.
if [ ! -w "${HOME:-/nonexistent}" ]; then
    export HOME="$PROHAP_WORK/.home"
    mkdir -p "$HOME"
fi
export XDG_CACHE_HOME="$PROHAP_WORK/.cache"
mkdir -p "$XDG_CACHE_HOME"

log "$image, work directory $PROHAP_WORK"
if [ "$command" = verify ]; then
    exec python3 "$PROHAP_HOME/verify/verify.py" "$@"
fi

if [ "$#" -eq 0 ]; then
    log "no arguments: showing a dry run. Run with e.g. '--cores 16 -p' to execute."
    set -- --dry-run
fi
log "snakemake $*"
# --rerun-triggers mtime: a newer image changes the workflow's code, and Snakemake would otherwise rerun every job
# whose outputs already exist. Before "$@", so that an explicit --rerun-triggers still wins.
exec snakemake \
    --snakefile "$PROHAP_HOME/Snakefile" \
    --directory "$PROHAP_WORK" \
    --rerun-triggers mtime \
    "$@"
