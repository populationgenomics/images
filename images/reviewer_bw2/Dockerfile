FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
ARG VERSION=${VERSION:-0.2.7-bw2-fork-0.2}
# Install Git to clone the repository
RUN apt-get update && \
    apt-get install -y \
    cmake \
    g++ \
    git \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    zlib1g-dev
# Clone the repository, jump into it
# build the REViewer executable, +x, move to /usr/bin
RUN git clone https://github.com/bw2/REViewer && \
    cd REViewer && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    chmod +x /REViewer/build/install/bin/REViewer && \
    mv /REViewer/build/install/bin/REViewer /usr/bin
# Set the entry point for the Docker container
ENTRYPOINT ["REViewer"]
