FROM ubuntu:20.04
MAINTAINER Alena Musilova, alena.musilova@lfmotol.cuni.cz

LABEL \
  version="v2.2" \
  description="stripy-pipeline"

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/London

# Run update and install necessary libraries, including BWA and python3
RUN apt-get update -y && apt-get install -y --no-install-recommends \
    build-essential \
    bwa \
    bzip2 \
    ca-certificates \
    gcc \
    git \
    libbz2-dev \
    libcurl4-openssl-dev \
    libffi-dev \
    liblzma-dev \
    libncurses5-dev \
    libssl-dev \
    make \
    python3 \
    python3-dev \
    python3-pip \
    software-properties-common \
    wget \
    zlib1g-dev

WORKDIR /usr/local/bin/

# Install Samtools
ENV SAMTOOLS_VERSION="1.16.1"
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure \
    && make \
    && make install \
    && make clean

WORKDIR /usr/local/bin/

# Install ExpansionHunter
ENV EXPANSIONHUNTER_VERSION="v5.0.0"
RUN wget https://github.com/Illumina/ExpansionHunter/releases/download/${EXPANSIONHUNTER_VERSION}/ExpansionHunter-${EXPANSIONHUNTER_VERSION}-linux_x86_64.tar.gz \
    && tar xvzf ExpansionHunter-${EXPANSIONHUNTER_VERSION}-linux_x86_64.tar.gz \
    && rm ExpansionHunter-${EXPANSIONHUNTER_VERSION}-linux_x86_64.tar.gz \
    && chmod 755 ExpansionHunter-${EXPANSIONHUNTER_VERSION}-linux_x86_64/bin/ExpansionHunter \
    && ln -s /usr/local/bin/ExpansionHunter-${EXPANSIONHUNTER_VERSION}-linux_x86_64/bin/ExpansionHunter /usr/bin/ExpansionHunter

# Install REViewer
ENV REVIEWER_VERSION="v0.2.7"
RUN wget https://github.com/Illumina/REViewer/releases/download/${REVIEWER_VERSION}/REViewer-${REVIEWER_VERSION}-linux_x86_64.gz \
    && gunzip REViewer-${REVIEWER_VERSION}-linux_x86_64.gz \
    && mv REViewer-${REVIEWER_VERSION}-linux_x86_64 /usr/local/bin/REViewer \
    && chmod 755 /usr/local/bin/REViewer \
    && ln -s /usr/local/bin/REViewer /usr/bin/REViewer

# Install STRipy by cloning from git repo (the most updated version)
RUN git clone https://gitlab.com/andreassh/stripy-pipeline.git \
    && chmod 755 stripy-pipeline/batch.sh \
    && python3 -m pip install -r stripy-pipeline/requirements.txt

WORKDIR /usr/local/bin/stripy-pipeline