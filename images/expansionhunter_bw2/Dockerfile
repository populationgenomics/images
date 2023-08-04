FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

ARG VERSION=${VERSION:-5.0.0-bw2}

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
    zlib1g-dev

# Clone the repository, jump into it
# build the ExpansionHunter executable, +x, move to /usr/bin
RUN git clone https://github.com/bw2/ExpansionHunter && \
    cd ExpansionHunter && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    chmod +x /ExpansionHunter/build/install/bin/ExpansionHunter && \
    mv /ExpansionHunter/build/install/bin/ExpansionHunter /usr/bin

# Set the entry point for the Docker container
ENTRYPOINT ["ExpansionHunter"]
