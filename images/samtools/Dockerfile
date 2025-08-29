FROM debian:bookworm-slim AS base

ENV VERSION=1.21

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
        libbz2-1.0 \
        libcurl4 \
        liblzma5 \
        libncurses5-dev \
        libssl3 \
        zlib1g && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

# second stage - add compile-relevant tools, obtain release, install from release
FROM base AS compiler
RUN apt-get update && apt-get install -y \
        bzip2 \
        gcc \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        make \
        wget \
        zlib1g-dev

	# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${VERSION}/samtools-${VERSION}.tar.bz2 && \
    tar -xf samtools-${VERSION}.tar.bz2 && \
    cd samtools-${VERSION} && \
    ./configure --enable-libcurl --enable-s3 --enable-gcs && \
    make && \
    make DESTDIR=/samtools_install install && \
    make -C htslib-$VERSION DESTDIR=/samtools_install install

FROM base
COPY --from=compiler /samtools_install/usr/local/bin/* /usr/local/bin/
