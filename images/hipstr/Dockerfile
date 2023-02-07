FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

# hadolint ignore=DL3008
RUN apt-get update -y && apt-get install -y --no-install-recommends -qqy \
    apt-transport-https \ 
    ca-certificates \
    git \
    libbz2-dev \
    liblzma-dev \
    make \
    libhts-dev \
    g++ \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev\
     && apt-get clean \
     && rm -rf /var/lib/apt/lists/*

RUN update-ca-certificates

# Install HipSTR
RUN git clone https://github.com/HipSTR-Tool/HipSTR && \
    cd HipSTR && \
    make

ENV PATH="/ExpansionHunter/bin/:$PATH"

# Clean up
RUN rm -rf /tmp/* \
    /var/tmp/* \
    /var/cache/apt/*