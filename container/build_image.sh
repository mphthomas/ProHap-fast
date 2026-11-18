#!/usr/bin/env bash
# Build the prohap-fast image of one Ensembl release, one platform at a time, with its tags and labels.
#
#   container/build_image.sh [--release 113] [--platform linux/amd64,linux/arm64] [--name prohap-fast]
#                            [--no-cache] [-- EXTRA DOCKER BUILDX BUILD ARGUMENTS]
#
# Tags, in the local image store:
#   <name>:<version>-ensembl.v<release>-<arch>   every platform built (arch: amd64, arm64)
#   <name>:<version>-ensembl.v<release>          the image of this machine's platform, when it is built
# <version> is the VERSION file. Labels: org.opencontainers.image.version (VERSION), .revision (the commit, with
# -dirty when the files the image copies have uncommitted changes) and .created (now, UTC).
#
# The references are built on this machine's platform, natively, and the image of the other platform reuses them
# from the build cache, so both images share the references layer. Build an image that is to be released with
# --no-cache, so that no layer comes from an older build; that rebuilds the references for each platform too (about
# 4 minutes each), and the two references layers then differ in the gzip timestamp of crap_tagged.fa.gz. Arguments
# after -- go to docker buildx build (CI passes --cache-from and --cache-to). Needs Docker with buildx; building
# for a platform other than this machine's needs emulation (Docker Desktop on Apple Silicon runs linux/amd64 with
# Rosetta).
set -euo pipefail

repo=$(cd "$(dirname "$0")/.." && pwd)
release=113
platforms=linux/amd64,linux/arm64
name=prohap-fast
no_cache=()
extra=()
while [ $# -gt 0 ]; do
    case "$1" in
        --release) release=$2; shift 2 ;;
        --platform) platforms=$2; shift 2 ;;
        --name) name=$2; shift 2 ;;
        --no-cache) no_cache=(--no-cache); shift ;;
        --) shift; extra=("$@"); break ;;
        -h | --help) sed -n '2,/^set -e/p' "$0" | sed '$d; s/^# \{0,1\}//'; exit 0 ;;
        *) echo "build_image.sh: unknown argument $1 (see --help)" >&2; exit 2 ;;
    esac
done
[[ $release =~ ^[0-9]+$ ]] || { echo "build_image.sh: --release must be a number, not '$release'" >&2; exit 2; }

version=$(tr -d '[:space:]' < "$repo/VERSION")
[[ $version =~ ^[0-9A-Za-z][0-9A-Za-z.-]*$ ]] || { echo "build_image.sh: VERSION '$version' cannot go in a tag" >&2; exit 2; }
tag="$name:$version-ensembl.v$release"

revision=unknown
if git -C "$repo" rev-parse --git-dir > /dev/null 2>&1; then
    revision=$(git -C "$repo" rev-parse HEAD)
    if [ -n "$(git -C "$repo" status --porcelain -- Dockerfile .dockerignore VERSION Snakefile entrypoint.sh crap.fasta \
               config_tre_template.yaml envs container legacy src tools verify)" ]; then
        revision="$revision-dirty"
    fi
fi
created=$(date -u +%Y-%m-%dT%H:%M:%SZ)
native="linux/$(docker version --format '{{.Server.Arch}}')"

IFS=, read -r -a platform_list <<< "$platforms"
for platform in "${platform_list[@]}"; do
    case "$platform" in
        linux/amd64 | linux/arm64) ;;
        *) echo "build_image.sh: unsupported platform $platform (linux/amd64, linux/arm64)" >&2; exit 2 ;;
    esac
    tags=(--tag "$tag-${platform#linux/}")
    [ "$platform" = "$native" ] && tags+=(--tag "$tag")
    echo "[build_image] $platform: ${tags[*]//--tag /}" >&2
    start=$(date +%s)
    docker buildx build "$repo" \
        --platform "$platform" \
        --build-arg ENSEMBL_RELEASE="$release" \
        --build-arg PROHAP_VERSION="$version" \
        --build-arg VCS_REF="$revision" \
        --build-arg BUILD_DATE="$created" \
        --provenance=false --sbom=false \
        "${tags[@]}" ${no_cache[@]+"${no_cache[@]}"} --load ${extra[@]+"${extra[@]}"}
    echo "[build_image] $platform built in $(( $(date +%s) - start )) s:" \
         "$(docker image inspect --format '{{.Size}}' "$tag-${platform#linux/}" | awk '{printf "%.2f GB", $1 / 1e9}')" >&2
done
