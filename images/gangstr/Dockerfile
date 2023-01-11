FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

ARG VERSION=${VERSION:-2.5.0}

# hadolint ignore=DL3008
RUN apt-get update && apt-get -y --no-install-recommends install wget && \
wget --progress=dot:giga --no-check-certificate https://github.com/gymreklab/GangSTR/releases/download/v2.5/GangSTR-${VERSION}-Source.tar.gz && \
gunzip GangSTR-${VERSION}-Source.tar.gz && \
chmod +x GangSTR-${VERSION}-Source.tar.gz && \
apt-get clean && \
rm -rf /var/lib/apt/lists/* \
    /tmp/* \
    /var/tmp/* \
    /var/cache/apt/* \
    GangSTR-${VERSION}-Source.tar.gz
