FROM australia-southeast1-docker.pkg.dev/cpg-common/images/nf_base:5.8.1-1 AS base

ARG VERSION=1.24

LABEL \
  version="v${VERSION}" \
  description="NF Base Samtools install"

RUN apt-get update && apt-get install -y --no-install-recommends \
        libbz2-1.0 \
        libcurl4 \
        liblzma5 \
        libncursesw6 \
        libssl3 \
        zlib1g && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

FROM base AS compiler

RUN apt-get update && apt-get install -y --no-install-recommends \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libncurses-dev \
        libssl-dev \
        zlib1g-dev && \
    wget https://github.com/samtools/samtools/releases/download/${VERSION}/samtools-${VERSION}.tar.bz2 && \
    tar -xf samtools-${VERSION}.tar.bz2 && \
    cd samtools-${VERSION} && \
    ./configure --enable-libcurl --enable-s3 --enable-gcs && \
    make && \
    strip samtools && \
    make DESTDIR=/samtools_install install && \
    make -C htslib-$VERSION DESTDIR=/samtools_install install

FROM base
COPY --from=compiler /samtools_install/usr/local/bin/* /usr/local/bin/
