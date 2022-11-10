FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
ARG VERSION=${VERSION:-5.0.0}

# hadolint ignore=DL3008
RUN apt-get update && apt-get install --no-install-recommends -qqy \
    python3-dev=3.8.2-0ubuntu2 \
    python3-pip=20.0.2-5ubuntu1.6 \
    python3.9\
    python-dev \
    build-essential \
    git \
    libbz2-dev \
    liblzma-dev \
    make \
    pkg-config \
    wget \
    unzip \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev\
     && apt-get clean \
     && rm -rf /var/lib/apt/lists/*

# Install ExpansionHunter
RUN wget --progress=dot:giga https://github.com/Illumina/ExpansionHunter/releases/download/${EH_VERSION}/ExpansionHunter-v${VERSION}-linux_x86_64.tar.gz && \
    tar xzf ExpansionHunter-v${VERSION}-linux_x86_64.tar.gz && \
    rm ExpansionHunter-v${VERSION}-linux_x86_64.tar.gz && \
    mv /ExpansionHunter-v${VERSION}-linux_x86_64 /ExpansionHunter
ENV PATH="/ExpansionHunter/bin/:$PATH"

# Clean up
RUN rm -rf /tmp/* \
    /var/tmp/* \
    /var/cache/apt/*
