ARG PY_VER=${PY_VER:-3.10}
ARG UV_VER=${UV_VER:-0.6.4}
ARG VERSION=${VERSION:-0.2.134.cpg1}

FROM ghcr.io/astral-sh/uv:${UV_VER}-python${PY_VER}-bookworm-slim AS basic

ENV HAIL_QUERY_BACKEND="service"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install --no-install-recommends -y \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        gnupg \
        openjdk-17-jdk-headless \
        software-properties-common \
        wget \
        zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    pip install --no-cache-dir --upgrade pip

FROM basic AS builder

ENV HAIL_QUERY_BACKEND="service"

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential \
        curl \
        g++ \
        gcc \
        git \
        libfontconfig \
        liblapack3 \
        libopenblas-dev \
        rsync && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

ENV UV_PROJECT_ENVIRONMENT=/usr/local

# Install Hail from the CPG fork as a wheel
RUN uv pip --no-cache-dir install --system \
        build && \
    # Install Hail from the CPG fork.
    git clone https://github.com/populationgenomics/hail.git && \
    cd hail && \
    git checkout $VERSION && \
    cd hail && \
    # Generate a Wheel, but don't install it in this layer
    # DEPLOY_REMOTE avoids a dev suffix being appended to dataproc initialization paths.
    make wheel DEPLOY_REMOTE=1 HAIL_RELEASE_MODE=1

FROM basic AS installer

# this RUN command binds in the wheel from the build layer without copying it in permanently
RUN --mount=type=bind,from=builder,source=/hail/hail/build/deploy/dist,target=/wheel \
    uv pip install --system --no-cache-dir \
        /wheel/hail-*-py3-none-any.whl
