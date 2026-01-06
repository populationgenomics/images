FROM debian:bookworm-slim

ENV MAMBA_ROOT_PREFIX /root/micromamba
ENV PATH $MAMBA_ROOT_PREFIX/bin:$PATH
ARG VERSION=2.4.6

RUN apt-get update && apt-get install -y git wget bash bzip2 zip && \
	wget -qO- https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba && \
    chmod +x /usr/local/bin/micromamba && \
    \
    mkdir -p ${MAMBA_ROOT_PREFIX} && \
    \
    # Perform the installation
    micromamba install -y --root-prefix ${MAMBA_ROOT_PREFIX} --prefix ${MAMBA_ROOT_PREFIX} \
    -c bioconda -c conda-forge \
    python=3.10 r-base=4.3.1 r-reshape2=1.4.4 r-tidyverse=2.0.0 \
    bioconductor-fraser=${VERSION} bioconductor-txdb.hsapiens.ucsc.hg38.knowngene bioconductor-org.hs.eg.db && \
    \
    # Clean up to keep the image small
    micromamba clean --all --yes && \
    rm -rf /var/lib/apt/lists/*
