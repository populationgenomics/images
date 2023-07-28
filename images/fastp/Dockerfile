FROM debian:bookworm-slim

ARG VERSION=${VERSION:-0.23.4}

RUN apt-get update && apt-get install -y git wget bash bzip2 zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    wget http://opengene.org/fastp/fastp.${VERSION} && \
    mv fastp.${VERSION} fastp && \
    chmod a+x ./fastp
