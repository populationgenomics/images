FROM ubuntu:22.04

ARG VERSION=${VERSION:-1.2.3}

LABEL \
  version="v${VERSION}" \
  description="mitoreport"

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/London

# Run update and install necessary libraries, mitoreport requires java11
RUN apt-get update -y && apt-get install -y --no-install-recommends \
    openjdk-11-jdk \
    apt-transport-https \
    build-essential \
    bzip2 \
    ca-certificates \
    curl \
    gcc \
    git \
    gnupg \
    gpg-agent \
    libbz2-dev \
    libcurl4-openssl-dev \
    libffi-dev \
    liblzma-dev \
    libncurses5-dev \
    libssl-dev \
    make \
    software-properties-common \
    wget \
    zlib1g-dev

## Install Google cloud-sdk
RUN curl https://sdk.cloud.google.com > install.sh && \
    bash install.sh --disable-prompts --install-dir=/opt && \
    rm install.sh

ENV PATH=$PATH:/opt/google-cloud-sdk/bin

WORKDIR /usr/local/bin/

# Instead of downloading JAR, clone and build
# Modify UI to always bypass file:// protocol check (to allow viewing from cloud storage)
RUN git clone https://github.com/bioinfomethods/mitoreport.git /tmp/mitoreport \
    && cd /tmp/mitoreport \
    && git checkout ${VERSION} \
    && sed -i "s/if (window\.location\.protocol === 'file:')/if (true || window.location.protocol === 'file:')/" ui/src/main.js \
    && ./gradlew build \
    && cp build/libs/mitoreport-*-all.jar /usr/local/bin/mitoreport.jar \
    && rm -rf /tmp/mitoreport

# Requires samtools for local cram > bam
ENV MAMBA_ROOT_PREFIX="/root/micromamba"
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH
RUN wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir ${MAMBA_ROOT_PREFIX} && \
    micromamba install -y --prefix ${MAMBA_ROOT_PREFIX} -c bioconda -c conda-forge \
    samtools && \
    rm -r /root/micromamba/pkgs
