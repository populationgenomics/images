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

# Install mitoreport and resourses
RUN curl -L -O https://github.com/bioinfomethods/mitoreport/releases/download/${VERSION}/mitoreport-${VERSION}-all.jar \
    && mv mitoreport-${VERSION}-all.jar mitoreport.jar \
    && curl -L -O https://github.com/bioinfomethods/mitoreport/releases/download/${VERSION}/resources.tgz \
    && tar -zxvf resources.tgz \
    && java -jar mitoreport.jar mito-map-download \
        --output resources/mito_map_annotations.json

# Requires samtools for local cram > bam
ENV MAMBA_ROOT_PREFIX="/root/micromamba"
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH
RUN wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir ${MAMBA_ROOT_PREFIX} && \
    micromamba install -y --prefix ${MAMBA_ROOT_PREFIX} -c bioconda -c conda-forge \
    samtools && \
    rm -r /root/micromamba/pkgs
