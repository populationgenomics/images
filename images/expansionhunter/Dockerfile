
FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
ENV PATH="/root/miniconda3/bin:${PATH}"
ENV EH_VERSION=v5.0.0

RUN apt-get update && apt-get install --no-install-recommends -qqy \
    python3-dev \
    python3-pip \
    python \
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
    libssl-dev

# Install ExpansionHunter
RUN wget https://github.com/Illumina/ExpansionHunter/releases/download/${EH_VERSION}/ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz && \
    tar xzf ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz && \
    rm ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz && \
    mv /ExpansionHunter-${EH_VERSION}-linux_x86_64 /ExpansionHunter
ENV PATH="/ExpansionHunter/bin/:$PATH"

# Clean up
RUN rm -rf /tmp/* \
    /var/tmp/* \
    /var/cache/apt/* \
    /var/lib/apt/lists/* \
    
