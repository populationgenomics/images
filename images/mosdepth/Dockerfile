FROM debian:bookworm-slim

ARG VERSION=${VERSION:-v0.3.6}

ENV PATH=$PATH:/opt/google-cloud-sdk/bin

WORKDIR /usr/local/bin/

RUN apt-get update && apt-get install -y git wget curl bash bzip2 zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    curl https://sdk.cloud.google.com > install.sh && \
    bash install.sh --disable-prompts --install-dir=/opt && \
    rm install.sh && \
    wget https://github.com/brentp/mosdepth/releases/download/${VERSION}/mosdepth && \
    wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig && \
    chmod a+x mosdepth wigToBigWig
