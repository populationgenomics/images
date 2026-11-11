FROM debian:bookworm-slim

ARG VERSION=2.10.0

# seqkit is a static binary. No USER is set, so the image runs as root: Hail Batch cannot run an image
# whose USER is a bare numeric UID (e.g. ontresearch/wf-common, USER 0).
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        procps \
        wget && \
    wget -q "https://github.com/shenwei356/seqkit/releases/download/v${VERSION}/seqkit_linux_amd64.tar.gz" \
        -O /tmp/seqkit.tar.gz && \
    tar -xzf /tmp/seqkit.tar.gz -C /usr/local/bin seqkit && \
    rm /tmp/seqkit.tar.gz && \
    apt-get purge -y wget && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*
