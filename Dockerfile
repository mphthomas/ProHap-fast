# syntax=docker/dockerfile:1.6
#
# Self-contained ProHap image for air-gapped Singularity TREs.
#
# Build on an internet-connected host (targeting linux/amd64 from a Mac):
#   docker buildx build --platform=linux/amd64 -t <repo>/prohap:<tag> --load .
#
# Override the Ensembl release at build time. This MUST match the
# ensembl_release / annotationFilename in the config.yaml you run with,
# otherwise the workflow will try to download the missing release and fail
# with no network:
#   docker buildx build --platform=linux/amd64 --build-arg ENSEMBL_RELEASE=114 ...

ARG ENSEMBL_RELEASE=113
ARG SNAKEMAKE_VERSION=7.32.4

# -------- Stage A: fetch Ensembl reference files (kept gzipped) --------
FROM --platform=$TARGETPLATFORM debian:bookworm-slim AS refs
ARG ENSEMBL_RELEASE
RUN apt-get update \
 && apt-get install -y --no-install-recommends wget ca-certificates \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /refs
RUN set -eux; \
    base="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}"; \
    wget -q "${base}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.chr_patch_hapl_scaff.gtf.gz"; \
    wget -q "${base}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz";  \
    wget -q "${base}/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz";    \
    wget -q "${base}/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz";    \
    # fail the build now rather than at run time inside the TRE
    for f in *.gz; do gzip -t "$f" || exit 1; done; \
    ls -l /refs

# -------- Stage B: runtime image --------
FROM --platform=$TARGETPLATFORM mambaorg/micromamba:1.5.10 AS runtime
ARG ENSEMBL_RELEASE
ARG SNAKEMAKE_VERSION

USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/* \
 && mkdir -p /opt/prohap /prohap/work \
 && chmod 0777 /prohap/work
USER $MAMBA_USER

# Runtime environment. Snakemake is pinned: v8+ reworked the executor and
# conda-deployment flags, and this workflow is written against v7 semantics.
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/prohap.yaml /tmp/env.yaml
RUN micromamba create -y -n ProHap -f /tmp/env.yaml \
 && micromamba install -y -n ProHap -c conda-forge -c bioconda \
      "snakemake-minimal=${SNAKEMAKE_VERSION}" pyyaml \
 && micromamba clean -a -y \
 && rm /tmp/env.yaml

USER root

COPY --from=refs /refs/ /opt/prohap/refs/

# Copy only what the workflow needs. Listing the paths explicitly keeps build
# artifacts, test output and editor state out of the image regardless of what
# happens to be lying around in the build context.
COPY Snakefile           /opt/prohap/Snakefile
COPY entrypoint.sh       /opt/prohap/entrypoint.sh
COPY src/                /opt/prohap/src/
COPY legacy/             /opt/prohap/legacy/
COPY envs/               /opt/prohap/envs/
COPY container/          /opt/prohap/container/
# crap.fasta only: contaminant sequences are generic and every run needs them.
# No example data, sample metadata or example configs are shipped - the image
# carries the tool and its references, nothing cohort-specific.
COPY crap.fasta          /opt/prohap/crap.fasta

# Singularity runs the container as the invoking host user, whose uid will not
# match any user in this image, so every path has to be world readable and the
# entrypoint world executable.
RUN chmod -R a+rX /opt/prohap /opt/conda \
 && chmod 0755 /opt/prohap/entrypoint.sh

USER $MAMBA_USER

# PROHAP_WORK is deliberately not set: the entrypoint defaults it to the
# current directory, which under Singularity is the caller's own bind-mounted
# directory. Set it explicitly only to override that.
#
# PYTHONHASHSEED: ProHap orders FASTA records by the Python hash of the
# sequence, and string hashing is randomised per process, so record order and
# the enshap_* accession numbering vary between otherwise identical runs. The
# content is unaffected, but pinning the seed makes a run byte-reproducible.
ENV ENSEMBL_RELEASE=${ENSEMBL_RELEASE} \
    PATH=/opt/conda/envs/ProHap/bin:$PATH \
    PROHAP_HOME=/opt/prohap \
    PROHAP_REFS=/opt/prohap/refs \
    PYTHONHASHSEED=0

LABEL org.opencontainers.image.title="ProHap" \
      org.opencontainers.image.description="ProHap/ProVar protein haplotype database builder, packaged for air-gapped TREs" \
      org.opencontainers.image.source="https://github.com/ProGenNo/ProHap" \
      io.prohap.ensembl_release="${ENSEMBL_RELEASE}" \
      io.prohap.snakemake_version="${SNAKEMAKE_VERSION}"

WORKDIR /prohap/work
ENTRYPOINT ["/opt/prohap/entrypoint.sh"]
