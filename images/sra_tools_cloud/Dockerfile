FROM debian:bookworm-slim AS installer

ENV VERSION=3.4.1

# The "-cloud" tarball ships a system-wide /etc/ncbi config with
# /config/default = "true", so the tools run for any user without the
# per-user GUID configuration step.
RUN apt-get update && apt-get install -y \
        ca-certificates \
        curl && \
    curl -fsSL "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${VERSION}/sratoolkit.${VERSION}-alma_linux64-cloud.tar.gz" -o sratoolkit.tar.gz && \
    mkdir /sratoolkit && \
    tar -xzf sratoolkit.tar.gz -C /sratoolkit && \
    printf '/repository/user/cache-disabled = "true"\n' >> /sratoolkit/etc/ncbi/default.kfg

FROM debian:bookworm-slim AS htslib_builder

ENV HTSLIB_VERSION=1.24

# Built from source because bookworm's tabix package ships htslib 1.16.
# --disable-libcurl: bgzip/tabix don't need remote file access, and it
# keeps libcurl out of the runtime dependencies.
RUN apt-get update && apt-get install -y \
        bzip2 \
        ca-certificates \
        curl \
        gcc \
        libbz2-dev \
        libdeflate-dev \
        liblzma-dev \
        make \
        zlib1g-dev && \
    curl -fsSL "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2" | \
        tar -xjf - && \
    cd "htslib-${HTSLIB_VERSION}" && \
    ./configure --disable-libcurl && \
    make -j"$(nproc)" && \
    make install

FROM debian:bookworm-slim

ENV VERSION=3.4.1

LABEL \
  version="v${VERSION}" \
  description="NCBI SRA Toolkit plus gcloud CLI, pigz, and bgzip, for fetching SRA objects from cloud buckets"

# libxml-libxml-perl is required by the perl-based loader tools;
# libbz2/liblzma/libdeflate are htslib runtime dependencies
RUN apt-get update && apt-get install -y \
        ca-certificates \
        curl \
        gnupg \
        libbz2-1.0 \
        libdeflate0 \
        liblzma5 \
        libxml-libxml-perl \
        pigz && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
        > /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
        gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    apt-get update && \
    apt-get install -y --no-install-recommends google-cloud-cli && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

COPY --from=installer /sratoolkit/usr/local/ncbi /usr/local/ncbi
COPY --from=installer /sratoolkit/etc/ncbi /etc/ncbi
COPY --from=htslib_builder /usr/local/bin/bgzip /usr/local/bin/tabix /usr/local/bin/
COPY --from=htslib_builder /usr/local/lib/libhts.so* /usr/local/lib/
RUN ldconfig

ENV PATH="/usr/local/ncbi/sra-tools/bin:${PATH}"
