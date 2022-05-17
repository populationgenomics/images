FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

# hadolint ignore=DL3008
WORKDIR /opt
RUN apt-get update && apt-get -y --no-install-recommends install wget && \
wget --progress=dot:giga --no-check-certificate https://github.com/Illumina/REViewer/releases/download/v0.2.7/REViewer-v0.2.7-linux_x86_64.gz && \
gunzip REViewer-v0.2.7-linux_x86_64.gz && \
chmod +x REViewer-v0.2.7-linux_x86_64 && \
ln REViewer-v0.2.7-linux_x86_64 /bin/reviewer && \
apt-get clean && \
rm -rf /var/lib/apt/lists/* \
    /tmp/* \
    /var/tmp/* \
    /var/cache/apt/* \
    REViewer-v0.2.7-linux_x86_64.gz
    