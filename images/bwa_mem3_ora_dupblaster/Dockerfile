FROM debian:trixie-slim AS builder

ARG VERSION=v0.15.0
ARG LLVM_VERSION=19
ARG ORAD_VERSION=2.7.0
ARG DUPBLASTER_VERSION=0.3.0
ARG BAZAM_VERSION=1.0.1

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update && apt-get install -y --no-install-recommends \
      autoconf \
      autoconf-archive \
      automake \
      build-essential \
      ca-certificates \
      clang-${LLVM_VERSION} \
      cmake \
      curl \
      git \
      libdeflate-dev \
      libomp-${LLVM_VERSION}-dev \
      libtool \
      pkg-config \
      zlib1g-dev && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

RUN git clone --depth 1 --branch ${VERSION} --recursive https://github.com/fg-labs/bwa-mem3.git /build && \
    cd /build && \
    make -j"$(nproc)" BASELINE_ARCH=avx2 CC=clang-${LLVM_VERSION} CXX=clang++-${LLVM_VERSION} && \
    install -Dm755 bwa-mem3 /out/bin/bwa-mem3

# Illumina ORA decompression
# https://sapac.support.illumina.com/downloads/dragen_decompression_software-linux.html
RUN curl -fsSL https://s3.amazonaws.com/webdata.illumina.com/downloads/software/dragen-decompression/orad.${ORAD_VERSION}.linux.tar.gz | \
        tar -xz -C /tmp && \
    install -Dm755 /tmp/orad.${ORAD_VERSION}.linux/orad /out/bin/orad

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | \
        sh -s -- -y --profile minimal && \
    /root/.cargo/bin/cargo install --locked --root /out dupblaster --version ${DUPBLASTER_VERSION}

# https://github.com/ssadedin/bazam#download-a-release-from-github
RUN curl -fsSL -o /out/share/bazam/bazam.jar --create-dirs \
        https://github.com/ssadedin/bazam/releases/download/${BAZAM_VERSION}/bazam.jar && \
    printf '#!/bin/sh\nexec java ${JAVA_OPTS} -jar /usr/local/share/bazam/bazam.jar "$@"\n' > /out/bin/bazam && \
    chmod 755 /out/bin/bazam

FROM debian:trixie-slim

ARG LLVM_VERSION=19

# procps provides ps, which Nextflow needs to record task CPU and memory use
RUN apt-get update && apt-get install -y --no-install-recommends \
      libdeflate0 \
      libomp5-${LLVM_VERSION} \
      openjdk-21-jre-headless \
      procps \
      zlib1g && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

COPY --from=builder /out/bin/bwa-mem3 /out/bin/orad /out/bin/dupblaster /out/bin/bazam /usr/local/bin/
COPY --from=builder /out/share/bazam/bazam.jar /usr/local/share/bazam/bazam.jar

RUN bwa-mem3 version && \
    orad --version && \
    dupblaster --version && \
    { bazam -h 2>&1 || true; } | grep -q 'usage: java -jar bazam.jar'
