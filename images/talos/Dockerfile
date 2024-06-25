FROM python:3.10-bullseye

ARG VERSION=${VERSION:-main}

RUN apt update && apt install -y \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        git \
        gnupg \
        openjdk-11-jdk-headless \
        wget \
        zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    # Google Cloud SDK: use the script-based installation, as the Debian package is outdated.
    curl https://sdk.cloud.google.com > install.sh && \
    bash install.sh --disable-prompts --install-dir=/opt && \
    rm install.sh

ENV PATH=$PATH:/opt/google-cloud-sdk/bin

RUN python3 -m pip install metamist git+https://github.com/populationgenomics/automated-interpretation-pipeline.git@${VERSION}
