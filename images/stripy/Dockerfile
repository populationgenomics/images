FROM

ARG VERSION=v2.5

LABEL \
  version="v${VERSION}" \
  description="stripy-pipeline"

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/London

# Run update and install necessary libraries, including BWA and python3
RUN apt-get update -y && apt-get install -y --no-install-recommends \
		apt-transport-https \
		build-essential \
		bwa \
		bzip2 \
		ca-certificates \
		curl \
		gcc \
		git \
		gnupg \
		gpg-agent \
		jq \
		libbz2-dev \
		libcurl4-openssl-dev \
		libffi-dev \
		liblzma-dev \
		libssl-dev \
		make \
		python3 \
		python3-dev \
		python3-pip \
		software-properties-common \
		wget \
		zlib1g-dev && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    pip install --no-cache-dir --upgrade pip

## Install Google cloud-sdk
RUN curl https://sdk.cloud.google.com > install.sh && \
    bash install.sh --disable-prompts --install-dir=/opt && \
    rm install.sh

ENV PATH=$PATH:/opt/google-cloud-sdk/bin

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
    && cd stripy-pipeline  \
    && git checkout $VERSION \
    && cd .. \
    && chmod 755 stripy-pipeline/batch.sh \
    && python3 -m pip install -r stripy-pipeline/requirements.txt

WORKDIR /usr/local/bin/stripy-pipeline

# Install gcsfuse to mount VEP cache files
RUN GCSFUSE_REPO=gcsfuse-$(lsb_release -c -s) && \
    echo "deb http://packages.cloud.google.com/apt $GCSFUSE_REPO main" | tee /etc/apt/sources.list.d/gcsfuse.list && \
    wget -qO - https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y gcsfuse
