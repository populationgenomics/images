FROM australia-southeast1-docker.pkg.dev/cpg-common/images-dev/nf_base:0.1.0-1 AS base

ARG VERSION=1.23.1

LABEL \
  version="v${VERSION}" \
  description="NF Base Bcftools install"

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
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libssl-dev \
        zlib1g-dev && \
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
