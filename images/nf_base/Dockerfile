FROM python:3.11-slim-trixie AS base

ARG VERSION=5.8.0

LABEL \
  version="v${VERSION}" \
  description="CPG NextFlow Base Image"

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Australia/Sydney

# Run update and install generic core packages
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        build-essential \
        bzip2 \
        ca-certificates \
        curl \
        git \
        gnupg \
        procps \
        wget && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    pip install --no-cache-dir \
    metamist \
    cpg-utils>=5.8.1
