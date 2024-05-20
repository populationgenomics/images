FROM debian:buster-slim

ENV MAMBA_ROOT_PREFIX /root/micromamba
ENV PATH $MAMBA_ROOT_PREFIX/bin:$PATH

ARG VERSION=${VERSION:-6.0.1}

# Install necessary packages
RUN apt-get update && apt-get install -y \
    git wget bash bzip2 zip && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/cache/apt/*

# Clone a specific branch from a GitHub repository
ARG GIT_REPO=https://github.com/populationgenomics/TRTools.git
ARG GIT_BRANCH=skip-if-monomorphic

RUN git clone -b ${GIT_BRANCH} ${GIT_REPO}

# Set the working directory
WORKDIR /TRTools
