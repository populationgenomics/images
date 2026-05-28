FROM python:3.11-slim-trixie AS base

ARG VERSION=0.1.0

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
    rm -r /var/cache/apt/*

RUN pip install \
    metamist \
    cpg-utils \
