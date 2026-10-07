FROM python:3.11-slim-bookworm AS basic

# ENV rather than ARG so the hail checkout below is immune to the CI
# `--build-arg VERSION=<tag>-<n>` override, which is not a valid git ref.
ENV VERSION=0.2.138.cpg2

ENV HAIL_QUERY_BACKEND="service"

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        curl \
        git \
        gnupg \
        jq \
        procps \
        rsync \
        software-properties-common \
        wget \
        zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    pip install --no-cache-dir --upgrade pip

# Hail requires Java 11, which Debian bookworm no longer packages, so install
# Eclipse Temurin 11 from the Adoptium apt repository.
# https://adoptium.net/installation/linux/
RUN wget -qO - https://packages.adoptium.net/artifactory/api/gpg/key/public | \
        gpg --dearmor -o /usr/share/keyrings/adoptium.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/adoptium.gpg] https://packages.adoptium.net/artifactory/deb bookworm main" \
        > /etc/apt/sources.list.d/adoptium.list && \
    apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
        temurin-11-jdk && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

FROM basic AS builder

ENV HAIL_QUERY_BACKEND="service"

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential \
        g++ \
        gcc \
        libfontconfig \
        liblapack3 \
        libopenblas-dev \
        openjdk-17-jre-headless && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

RUN pip --no-cache-dir install \
        build && \
    git clone https://github.com/populationgenomics/hail.git && \
    cd hail && \
    git checkout $VERSION && \
    cd hail && \
    # Generate a Wheel, but don't install it in this layer
    # DEPLOY_REMOTE avoids a dev suffix being appended to dataproc initialization paths.
    make wheel DEPLOY_REMOTE=1 HAIL_BUILD_MODE=Release

FROM basic AS installer

# wheel needs to have a descriptive name, e.g. hail-0.2.134-py3-none-any.whl
# but this isn't knowable exactly from the hash/branch, so we use a wildcard
# this RUN command binds in the wheel from the build layer without copying it in permanently
RUN --mount=type=bind,from=builder,source=/hail/hail/build/deploy/dist,target=/wheel \
    pip install --no-cache-dir \
        /wheel/hail-*-py3-none-any.whl \
