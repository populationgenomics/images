FROM debian:bookworm-slim AS base

ENV VERSION=1.23

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
        libbz2-1.0 \
        libcurl4 \
        liblzma5 \
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
        libssl-dev \
        make \
        wget \
        zlib1g-dev && \
    # Install bcftools
    wget https://github.com/samtools/bcftools/releases/download/${VERSION}/bcftools-${VERSION}.tar.bz2 && \
    tar -xf bcftools-${VERSION}.tar.bz2 && \
    cd bcftools-${VERSION} && \
    # Remove any existing array to vcf plugins and replace with the latest from the gtc2vcf repository
    /bin/rm -f plugins/{idat2gtc.c,gtc2vcf.{c,h},affy2vcf.c} && \
    wget -P plugins \
    https://raw.githubusercontent.com/freeseek/gtc2vcf/master/idat2gtc.c \
    https://raw.githubusercontent.com/freeseek/gtc2vcf/master/gtc2vcf.c \
    https://raw.githubusercontent.com/freeseek/gtc2vcf/master/gtc2vcf.h \
    https://raw.githubusercontent.com/freeseek/gtc2vcf/master/affy2vcf.c \
    https://raw.githubusercontent.com/freeseek/gtc2vcf/master/BAFregress.c && \
    ./configure --enable-libcurl --enable-s3 --enable-gcs && \
    make && \
    strip bcftools plugins/*.so && \
    make DESTDIR=/bcftools_install install && \
    make -C htslib-$VERSION DESTDIR=/bcftools_install install


FROM base
COPY --from=compiler /bcftools_install/usr/local/bin/* /usr/local/bin/
COPY --from=compiler /bcftools_install/usr/local/libexec/bcftools/* /usr/local/libexec/bcftools/
