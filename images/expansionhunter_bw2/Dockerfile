FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

ARG VERSION=${VERSION:-b85cba004fba4e0223e632fae576bc7fad0da804}

# Install Git to clone the repository
RUN apt-get update && \
    apt-get install -y \
    cmake \
    g++ \
    git \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libssl-dev \
    zlib1g-dev \
    && \
    apt-get clean -y && \
    rm -rf /var/lib/apt/lists/*

# Clone the repository, jump into it
# build the ExpansionHunter executable, +x, move to /usr/bin
RUN git clone https://github.com/bw2/ExpansionHunter && \
    cd ExpansionHunter && \
    git checkout $VERSION && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    chmod +x /ExpansionHunter/build/install/bin/ExpansionHunter && \
    mv /ExpansionHunter/build/install/bin/ExpansionHunter /usr/bin && \
    rm -rf /ExpansionHunter

