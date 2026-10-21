FROM debian:bookworm-slim

ARG VERSION=0.15.3

RUN apt-get update && apt-get install -y --no-install-recommends \
      bowtie2 \
      ca-certificates \
      libgd-perl \
      procps \
      wget \
    && \
    rm -rf /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

RUN wget -q -O /tmp/fastq_screen.tar.gz \
    https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v${VERSION}.tar.gz && \
    tar -xzf /tmp/fastq_screen.tar.gz -C /opt/ && \
    ln -s /opt/FastQ-Screen-${VERSION}/fastq_screen /usr/local/bin/fastq_screen && \
    rm /tmp/fastq_screen.tar.gz
