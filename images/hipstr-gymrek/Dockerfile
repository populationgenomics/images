FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
ARG VERSION=${VERSION:-1.0.0}

# hadolint ignore=DL3008
RUN apt-get update -y && apt-get install -y --no-install-recommends -qqy \
    apt-transport-https \ 
    ca-certificates \
    g++ \
    git \
    libbz2-dev \
    libcurl4-openssl-dev \
    libhts-dev \
    liblzma-dev \
    libssl-dev\
    make \
    wget \
    zlib1g-dev \
     && apt-get clean \
     && rm -rf /var/lib/apt/lists/* \
     && rm -rf /tmp/* \
        /var/tmp/* \
        /var/cache/apt/*

RUN update-ca-certificates

# Install HipSTR, copy, make executable delete
RUN wget https://github.com/gymrek-lab/HipSTR/archive/refs/tags/v${VERSION}.tar.gz  && \
    tar xzf v${VERSION}.tar.gz && \ 
    rm v${VERSION}.tar.gz && \ 
    mv /HipSTR-1.0.0 /HipSTR 
    #cp HipSTR /usr/local/bin/ && \
    #chmod a+x HipSTR /usr/local/bin/ && \
    #cd .. && \
    #rm -rf HipSTR