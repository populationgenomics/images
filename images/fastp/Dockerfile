FROM debian:bookworm-slim

ARG VERSION=${VERSION:-0.23.4}

RUN apt-get update && apt-get install -y git wget bash bzip2 zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    wget http://opengene.org/fastp/fastp.${VERSION} && \
    chmod a+x fastp.${VERSION} && \
    ln -s /fastp.${VERSION} /usr/local/bin/fastp
