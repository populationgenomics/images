FROM debian:bookworm-slim AS base

ARG VERSION=${VERSION:-0.3.11}

SHELL ["/bin/bash", "-eo", "pipefail", "-c"]

RUN apt-get update && apt-get install -y \
        curl && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    curl -L -o '/usr/local/bin/mosdepth' https://github.com/brentp/mosdepth/releases/download/v${VERSION}/mosdepth && \
    chmod a+x /usr/local/bin/mosdepth
