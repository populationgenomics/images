FROM debian:trixie-slim AS base

ENV DEBIAN_FRONTEND=noninteractive
ARG VERSION=${VERSION:-0.3.1}

# Install dependencies
RUN apt-get update && apt-get install -y \
    git wget bzip2 zip curl bc && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*
# Install somalier
RUN wget https://github.com/brentp/somalier/releases/download/v${VERSION}/somalier && \
    mv ./somalier /usr/local/bin/somalier && \
    chmod +x /usr/local/bin/somalier
