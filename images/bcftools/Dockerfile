FROM debian:bookworm-slim as base

ARG VERSION=${VERSION:-1.20}

RUN apt-get update && apt-get install -y \
        bzip2 \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        zlib1g-dev && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

# second stage - add compile-relevant tools, obtain release, install from release
FROM base as compiler

RUN apt-get update && apt-get install -y \
        gcc \
        make \
        wget && \
    wget https://github.com/samtools/bcftools/releases/download/${VERSION}/bcftools-${VERSION}.tar.bz2 && \
    tar -xf bcftools-${VERSION}.tar.bz2 && \
    cd bcftools-${VERSION} && \
    ./configure --prefix=/bcftools_install && \
    make install

# the actual build - on minimal install, copy content from bcftools install directory
from base

COPY --from=compiler /bcftools_install/bin/* /usr/bin
