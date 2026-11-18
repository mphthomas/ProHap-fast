# syntax=docker/dockerfile:1.17.1@sha256:38387523653efa0039f8e1c89bb74a30504e76ee9f565e25c9a09841f9427b05
#
# prohap-fast release image: the workflow, its locked environment and every reference derived from one Ensembl
# release, pre-built. One image per release, tagged prohap-fast:<version>-ensembl.v<release>. The Ensembl files are
# downloaded while the image is built; nothing is downloaded when it runs.
#
# Build with container/build_image.sh, which sets the tags and the labels:
#
#   container/build_image.sh                                    # release 113, linux/amd64 and linux/arm64
#   container/build_image.sh --release 114 --platform linux/amd64 --no-cache
#
# Stages:
#   download   build platform    the four Ensembl downloads of the release, checked with gzip -t
#   env-build  build platform    the locked environment, so that the references are built natively
#   prebuild   build platform    container/prebuild_refs.py with legacy/src -> /refs/ensembl.v<release>/
#   env        target platform   the same locked environment
#   runtime    target platform   env + references + workflow
# The references are built once, on the build platform, and the same layer goes into the image of every platform.
# env-build and env run the same commands, so they are one set of layers when the two platforms are the same.

ARG ENSEMBL_RELEASE=113
ARG SNAKEMAKE_VERSION=7.32.4

# -------- download: the Ensembl files of the release, kept gzipped --------
FROM --platform=$BUILDPLATFORM debian:bookworm-slim@sha256:734a72b3d6142a71993b4bc580ca01d8ca687f77a2c422d96a961d3376b179bc AS download
ARG ENSEMBL_RELEASE
RUN apt-get update \
 && apt-get install -y --no-install-recommends wget ca-certificates \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /downloads
RUN set -eux; \
    base="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}"; \
    wget -q "${base}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.chr_patch_hapl_scaff.gtf.gz"; \
    wget -q "${base}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"; \
    wget -q "${base}/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz"; \
    wget -q "${base}/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"; \
    for f in *.gz; do gzip -t "$f"; done; \
    ls -l

# -------- env-build: the locked environment on the build platform --------
# Keep this stage and `env` identical apart from the FROM line.
FROM --platform=$BUILDPLATFORM mambaorg/micromamba:1.5.10@sha256:e3797091302382ea841498bc93a7b0a50f7c1448333d5e946d2d1608d0c5f43d AS env-build
ARG SNAKEMAKE_VERSION
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/prohap.linux-64.lock envs/prohap.linux-aarch64.lock envs/prohap.pip.lock \
     container/install_env.sh /tmp/envs/
RUN bash /tmp/envs/install_env.sh /tmp/envs "$SNAKEMAKE_VERSION" && rm -rf /tmp/envs
# pandas and pyarrow are pip wheels whose C++ extensions find libstdc++ through the system search path. On
# linux/arm64 that loads the base image's older copy first, and the environment's ICU (loaded by sqlite3, so by
# gffutils) then fails to import. On linux/amd64 numpy already loads the environment's copy, so this changes nothing
# there. The entrypoint sets both again, because Singularity does not reliably pass ENV through.
ENV PATH=/opt/conda/envs/ProHap/bin:$PATH \
    LD_LIBRARY_PATH=/opt/conda/envs/ProHap/lib

# -------- prebuild: every reference the workflow reads, built by the frozen scripts --------
# A stage built FROM another stage takes its files from that stage but its platform from its own FROM line, which
# defaults to the target platform. Without --platform here, the build of each target platform would run this stage
# (natively, but as a different cache entry) and the images would not share the references layer.
FROM --platform=$BUILDPLATFORM env-build AS prebuild
ARG ENSEMBL_RELEASE
USER root
COPY --from=download /downloads/ /downloads/
COPY legacy/src/ /build/legacy/src/
COPY container/prebuild_refs.py /build/container/prebuild_refs.py
COPY crap.fasta /build/crap.fasta
# --src and --contaminants default to /build/legacy/src and /build/crap.fasta. The result is world-readable.
RUN python3 /build/container/prebuild_refs.py --release "$ENSEMBL_RELEASE" --downloads /downloads --out /refs

# -------- env: the locked environment on the target platform --------
FROM mambaorg/micromamba:1.5.10@sha256:e3797091302382ea841498bc93a7b0a50f7c1448333d5e946d2d1608d0c5f43d AS env
ARG SNAKEMAKE_VERSION
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/prohap.linux-64.lock envs/prohap.linux-aarch64.lock envs/prohap.pip.lock \
     container/install_env.sh /tmp/envs/
RUN bash /tmp/envs/install_env.sh /tmp/envs "$SNAKEMAKE_VERSION" && rm -rf /tmp/envs
ENV PATH=/opt/conda/envs/ProHap/bin:$PATH \
    LD_LIBRARY_PATH=/opt/conda/envs/ProHap/lib

# -------- runtime --------
FROM env AS runtime
ARG ENSEMBL_RELEASE
ARG SNAKEMAKE_VERSION
USER root
# /prohap/work is the default working directory for runs that do not set one (docker run without -w).
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/* \
 && mkdir -p /opt/prohap/refs /prohap/work \
 && chmod 0777 /prohap/work

# The image runs as whichever user calls it (docker run -u, Singularity), so everything is world-readable. Modes are
# set as the files are copied: a chmod -R afterwards would store every file a second time in a new layer.
COPY --from=prebuild --chmod=a+rX /refs/ensembl.v${ENSEMBL_RELEASE}/ /opt/prohap/refs/ensembl.v${ENSEMBL_RELEASE}/
# Only what the workflow needs, path by path, so that nothing lying around in the build context gets in. No example
# data or sample metadata: the image carries the tool and its references, nothing cohort-specific. crap.fasta is the
# default contaminants_fasta; config_tre_template.yaml is the config to start from.
COPY --chmod=a+rX envs/ /opt/prohap/envs/
COPY --chmod=a+rX legacy/ /opt/prohap/legacy/
COPY --chmod=a+rX tools/ /opt/prohap/tools/
COPY --chmod=a+rX container/ /opt/prohap/container/
COPY --chmod=a+rX crap.fasta config_tre_template.yaml VERSION /opt/prohap/
COPY --chmod=a+rX verify/ /opt/prohap/verify/
COPY --chmod=a+rX entrypoint.sh Snakefile /opt/prohap/
COPY --chmod=a+rX src/ /opt/prohap/src/

# ENSEMBL_RELEASE names the release of the references in this image; the entrypoint refuses a config for another
# one. PYTHONHASHSEED: ProHap orders FASTA records by the Python hash of the sequence, so the seed must be fixed for
# byte-identical output (the Snakefile also sets it for every rule). PYTHONNOUSERSITE: Singularity mounts $HOME, and a
# ~/.local site-packages there must not shadow the environment's packages, in the rules (the Snakefile sets it too)
# or in Snakemake and verify themselves. The entrypoint sets all of these again, because Singularity passes the
# host's environment in and does not reliably apply the image's.
ENV ENSEMBL_RELEASE=${ENSEMBL_RELEASE} \
    PROHAP_HOME=/opt/prohap \
    PROHAP_REFS=/opt/prohap/refs \
    PYTHONHASHSEED=0 \
    PYTHONNOUSERSITE=1

USER $MAMBA_USER
# Fail the build, not a run: the environment must load what the workflow imports, and another user must be able to
# read every file of the environment and the workflow, and open the references read-only.
RUN set -eu; \
    python -c "import pandas, sqlite3, gffutils, Bio, vcf, yaml, pyarrow"; \
    bcftools --version > /dev/null; \
    unreadable="$(find /opt/prohap /opt/conda/envs/ProHap \( -type d ! -perm -o=rx \) -o \( ! -type d ! -type l ! -perm -o=r \) | head)"; \
    if [ -n "$unreadable" ]; then echo "not world-readable:" >&2; echo "$unreadable" >&2; exit 1; fi; \
    refs="/opt/prohap/refs/ensembl.v${ENSEMBL_RELEASE}"; \
    test "$(cat "$refs/release.txt")" = "$ENSEMBL_RELEASE"; \
    python -c "import sqlite3, sys; con = sqlite3.connect('file:' + sys.argv[1] + '?mode=ro', uri=True); \
assert con.execute('select count(*) from features').fetchone()[0] > 0" \
        "$refs/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.chr_patch_hapl_scaff_chr22.db"

# Labels last, so that a new build date or commit does not invalidate any layer. build_image.sh sets the arguments.
ARG PROHAP_VERSION=unknown
ARG VCS_REF=unknown
ARG BUILD_DATE=unknown
LABEL org.opencontainers.image.title="prohap-fast" \
      org.opencontainers.image.description="ProHap protein haplotype databases for biobank cohorts, with the references of Ensembl release ${ENSEMBL_RELEASE} pre-built" \
      org.opencontainers.image.source="https://github.com/mphthomas/ProHap-fast" \
      org.opencontainers.image.url="https://github.com/mphthomas/ProHap-fast" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.version="${PROHAP_VERSION}" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.created="${BUILD_DATE}" \
      io.prohap.ensembl_release="${ENSEMBL_RELEASE}" \
      io.prohap.snakemake_version="${SNAKEMAKE_VERSION}"

WORKDIR /prohap/work
ENTRYPOINT ["/opt/prohap/entrypoint.sh"]
