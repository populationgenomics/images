FROM debian:bullseye-slim

ENV MAMBA_ROOT_PREFIX /root/micromamba
ENV PATH $MAMBA_ROOT_PREFIX/bin:$PATH
ARG VERSION=${VERSION:-1.3.0}

RUN apt-get update && apt-get install -y \
        git wget bash bzip2 zip build-essential libboost-all-dev libz-dev && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir ${MAMBA_ROOT_PREFIX} && \
    micromamba install -y --prefix ${MAMBA_ROOT_PREFIX} -c bioconda -c conda-forge \
        bazam biobambam google-cloud-sdk samtools \
    && rm -r /root/micromamba/pkgs && \
    git clone --depth 1 --branch ${VERSION} https://github.com/populationgenomics/DRAGMAP && \
    HAS_GTEST=0 make -C DRAGMAP -j 8 install && \
    rm -r DRAGMAP
