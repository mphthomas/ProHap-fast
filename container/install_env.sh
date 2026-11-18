#!/usr/bin/env bash
# Create the ProHap environment from the locks, for the architecture this runs on. Used by the Dockerfile:
#
#   install_env.sh LOCK_DIR SNAKEMAKE_VERSION
#
# LOCK_DIR holds prohap.linux-64.lock, prohap.linux-aarch64.lock and prohap.pip.lock (envs/). The architecture comes
# from `uname -m` rather than a build argument, so the Dockerfile's two environment stages (one for the build
# platform, one for the target) run the same command and share their cache when the two platforms are the same.
set -euo pipefail
locks=$1
snakemake_version=$2

case "$(uname -m)" in
    x86_64) subdir=linux-64 ;;
    aarch64 | arm64) subdir=linux-aarch64 ;;
    *) echo "install_env.sh: no lock for architecture $(uname -m)" >&2; exit 1 ;;
esac

micromamba create -y -n ProHap -f "$locks/prohap.$subdir.lock"
# pandas and pyarrow come as wheels, hash-checked; nothing else is resolved.
micromamba run -n ProHap python -m pip install --no-deps --only-binary=:all: --no-cache-dir \
    --require-hashes -r "$locks/prohap.pip.lock"
# Snakemake is pinned in the locks: v8+ reworked the executor and conda-deployment flags, and this workflow is
# written against v7. The check keeps the image's io.prohap.snakemake_version label honest.
test "$(micromamba run -n ProHap snakemake --version)" = "$snakemake_version"
micromamba clean -a -y
