FROM debian:buster-slim

ARG VERSION=${VERSION:-1.20}


RUN apt-get update && apt-get install -y \
        bash \
        bzip2 \
        gcc \
        git \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        make \
        wget \
        zip \
        zlib1g-dev && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    git clone --recurse-submodules https://github.com/samtools/htslib.git --branch ${VERSION} --depth 1 && \
    git clone https://github.com/samtools/bcftools.git --branch ${VERSION} --depth 1 && \
    cd bcftools && \
    make && \
    mv bcftools /usr/local/bin && \
    cd - && \
    rm -rf bcftools htslib
