FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

ARG VERSION=${VERSION:-0.2.7}

# hadolint ignore=DL3008
RUN apt-get update && apt-get -y --no-install-recommends install wget && \
wget --progress=dot:giga --no-check-certificate https://github.com/Illumina/REViewer/releases/download/v${VERSION}/REViewer-v${VERSION}-linux_x86_64.gz && \
gunzip REViewer-v${VERSION}-linux_x86_64.gz && \
chmod +x REViewer-v${VERSION}-linux_x86_64 && \
apt-get clean && \
rm -rf /var/lib/apt/lists/* \
    /tmp/* \
    /var/tmp/* \
    /var/cache/apt/* \
    REViewer-v${VERSION}-linux_x86_64.gz
