FROM docker.io/docker/hardened-debian:12-slim AS base

ARG VERSION=${VERSION:-1.21}
ENV VERSION=1.21
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
        libbz2-1.0 \
        libcurl4 \
        liblzma5 \
        libssl3 \
        zlib1g && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

FROM base AS compiler
RUN apt-get update && apt-get install -y \
        bzip2 \
        gcc \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libssl-dev \
        make \
        wget \
        zlib1g-dev && \
    # Install htslib first
    wget https://github.com/samtools/htslib/releases/download/${VERSION}/htslib-${VERSION}.tar.bz2 && \
    tar -xf htslib-${VERSION}.tar.bz2 && \
    cd htslib-${VERSION} && \
    ./configure --enable-libcurl --enable-s3 --enable-gcs && \
    make && \
    make DESTDIR=/htslib_install install && \
    cd .. && \
    # Then install bcftools
    wget https://github.com/samtools/bcftools/releases/download/${VERSION}/bcftools-${VERSION}.tar.bz2 && \
    tar -xf bcftools-${VERSION}.tar.bz2 && \
    cd bcftools-${VERSION} && \
    ./configure --enable-libcurl --enable-s3 --enable-gcs && \
    make && \
    strip bcftools plugins/*.so && \
    make DESTDIR=/bcftools_install install

FROM base
COPY --from=compiler /htslib_install/usr/local/bin/* /usr/local/bin/
COPY --from=compiler /bcftools_install/usr/local/bin/* /usr/local/bin/
COPY --from=compiler /bcftools_install/usr/local/libexec/bcftools/* /usr/local/libexec/bcftools/
