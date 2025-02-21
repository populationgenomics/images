FROM debian:bookworm-slim AS base

#ARG VERSION=${VERSION:-1.04.00}

SHELL ["/bin/bash", "-eo", "pipefail", "-c"]

 RUN apt-get update && apt-get install -y curl && \
     rm -r /var/lib/apt/lists/* && \
     rm -r /var/cache/apt/* && \
     curl -L -o '/usr/local/bin/bigBedToBed' http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed && \
     chmod a+x /usr/local/bin/bigBedToBed