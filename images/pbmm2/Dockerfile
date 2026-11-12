FROM debian:bookworm-slim

ARG VERSION=26.2.0
ARG MAMBA_ROOT_PREFIX=/opt/micromamba

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# procps provides ps, which Nextflow needs to record task CPU and memory use
RUN apt-get update && apt-get install -y --no-install-recommends \
        bzip2 \
        ca-certificates \
        procps \
        wget && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/2.9.0/linux-64/micromamba-2.9.0-0.tar.bz2 | \
        tar -xj -C /usr/local bin/micromamba && \
    micromamba create -y --prefix /opt/pbmm2 -c conda-forge -c bioconda pbmm2=${VERSION} && \
    micromamba clean -y --all

ENV PATH=/opt/pbmm2/bin:$PATH
