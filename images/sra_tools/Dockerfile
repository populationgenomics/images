FROM debian:bookworm-slim AS installer

ARG VERSION=3.4.1

RUN apt-get update && apt-get install -y \
        ca-certificates \
        curl && \
    curl -s "https://ftp.ncbi.nlm.nih.gov/sra/sdk/${VERSION}/sratoolkit.${VERSION}-ubuntu64.tar.gz" -o sratoolkit.tar.gz && \
    tar -xzf sratoolkit.tar.gz && \
    mv "sratoolkit.${VERSION}-ubuntu64" /usr/local/sratoolkit

FROM debian:bookworm-slim

ARG VERSION

LABEL \
  version="v${VERSION}" \
  description="NCBI SRA Toolkit (prefetch, fasterq-dump, etc.)"

RUN apt-get update && apt-get install -y \
        ca-certificates && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

COPY --from=installer /usr/local/sratoolkit /usr/local/sratoolkit

ENV PATH="/usr/local/sratoolkit/bin:${PATH}"

# The tools refuse to run until a configuration with a GUID exists; also keep
# the download cache out of the image layer's home directory.
RUN mkdir -p /root/.ncbi && \
    printf '/LIBS/GUID = "%s"\n/repository/user/cache-disabled = "true"\n' \
        "$(cat /proc/sys/kernel/random/uuid)" > /root/.ncbi/user-settings.mkfg && \
    prefetch --version && \
    fasterq-dump --version
