FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

# hadolint ignore=DL3008
RUN apt update && apt -y install wget && \
wget https://github.com/Illumina/REViewer/releases/download/v0.2.7/REViewer-v0.2.7-linux_x86_64.gz && \
gunzip REViewer-v0.2.7-linux_x86_64.gz && \
chmod +x REViewer-v0.2.7-linux_x86_64 && \
apt-get clean && \
rm -rf /var/lib/apt/lists/*

# Clean up
RUN rm -rf /tmp/* \
    /var/tmp/* \
    /var/cache/apt/*


