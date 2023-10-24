FROM debian:bookworm-slim

ARG VERSION=${VERSION:-1.0.1}

RUN apt-get update && apt-get install -y git wget bash bzip2 zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    mkdir -p /tmp/sambamba && \
    cd /tmp/sambamba && \
    wget https://github.com/biod/sambamba/releases/download/v${VERSION}/sambamba-${VERSION}-linux-amd64-static.gz && \
    gunzip sambamba-${VERSION}-linux-amd64-static.gz && \
    mv sambamba-${VERSION}-linux-amd64-static /usr/local/bin/sambamba && \
    chmod a+x /usr/local/bin/sambamba
