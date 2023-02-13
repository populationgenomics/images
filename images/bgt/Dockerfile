FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
ARG VERSION=${VERSION:-v1.0}

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
    zlib1g-dev \
     && apt-get clean \
     && rm -rf /var/lib/apt/lists/* \
     && rm -rf /tmp/* \
        /var/tmp/* \
        /var/cache/apt/* \
    && update-ca-certificates \
    && git clone -b ${VERSION} https://github.com/lh3/bgt.git \
    && cd bgt \
    && make \
    && cp bgt /usr/local/bin/ \
    && chmod a+x /usr/local/bin/bgt \
    && cd .. \
    && rm -rf bgt