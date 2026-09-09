FROM debian:bookworm-slim

ARG VERSION=0.15.3

RUN apt-get update && apt-get install -y --no-install-recommends \
    perl wget ca-certificates bowtie2 libgd-perl && \
    rm -rf /var/lib/apt/lists/*

RUN wget -q -O /tmp/fastq_screen.tar.gz \
    https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v${VERSION}.tar.gz && \
    tar -xzf /tmp/fastq_screen.tar.gz -C /opt/ && \
    ln -s /opt/FastQ-Screen-${VERSION}/fastq_screen /usr/local/bin/fastq_screen && \
    rm /tmp/fastq_screen.tar.gz