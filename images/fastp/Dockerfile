FROM debian:bookworm-slim

ARG VERSION=${VERSION:-1.0.1}

RUN apt-get update && apt-get install -y \
      bzip2 \
      procps \
      wget \
      zip \
    && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    wget http://opengene.org/fastp/fastp.${VERSION} && \
    chmod a+x fastp.${VERSION} && \
    ln -s /fastp.${VERSION} /usr/local/bin/fastp
