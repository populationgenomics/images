FROM python:3.10-slim

ARG VERSION=${VERSION:-v.1.10}

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        autoconf \
        bedtools \
        build-essential \
        curl \
        git \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libncurses5-dev \
        python3-dev \
        tabix \
        wget && \
    apt-get clean -y  && \
    rm -rf /var/lib/apt/lists/* && \
    git clone https://github.com/andrewSharo/StrVCTVRE && \
    cd StrVCTVRE && \
    git checkout ${VERSION} && \
    pip install \
    cyvcf2 \
    Cython \
    joblib \
    numpy \
    pandas \
    pybedtools \
    pybigwig \
    pysam \
    scikit-learn==1.0.2

WORKDIR /StrVCTVRE
