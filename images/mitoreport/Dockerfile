FROM ubuntu:22.04

ARG VERSION=${VERSION:-1.0.0-beta-1}

LABEL \
  version="v${VERSION}" \
  description="mitoreport"

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/London

# Run update and install necessary libraries, including BWA and python3
RUN apt-get update -y && apt-get install -y --no-install-recommends \
    default-jre \
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

# Install mitoreport and resourses
RUN wget https://github.com/bioinfomethods/mitoreport/releases/download/${VERSION}/mitoreport-${VERSION}-all.jar  \
    && wget https://bioinfomethods.github.io/mitoreport/resources/resources.tgz \
    && tar -zxvf resources.tgz \
    && rm resources.tgz

ENV MAMBA_ROOT_PREFIX /root/micromamba
ENV PATH $MAMBA_ROOT_PREFIX/bin:$PATH

RUN wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir ${MAMBA_ROOT_PREFIX} && \
    micromamba install -y --prefix ${MAMBA_ROOT_PREFIX} -c bioconda -c conda-forge \
    samtools && \
    rm -r /root/micromamba/pkgs
