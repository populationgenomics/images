ARG VERSION=4.6.2.0

FROM australia-southeast1-docker.pkg.dev/cpg-common/images/gatk:4.6.2.0-21-ge8c49f600-NIGHTLY-SNAPSHOT AS base

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        curl \
        git \
        gnupg \
        jq \
        libgomp1 \
        procps \
        rsync \
        software-properties-common \
        unzip \
        wget \
        zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

# build bcftools from source
FROM base AS compiler

ARG BCFTOOLS_VERSION=1.22
RUN apt-get update && apt-get install -y \
        bzip2 \
        gcc \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libssl-dev \
        make \
        wget \
        zlib1g-dev && \
    # Install bcftools
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -xf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    # Reset PATH to exclude Conda, and unset the include variables
    export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin && \
    unset CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH LIBRARY_PATH && \
    # Run configure with explicit help from pkg-config
    PKG_CONFIG_PATH=/usr/lib/x86_64-linux-gnu/pkgconfig ./configure --enable-libcurl --enable-s3 --enable-gcs && \
    make && \
    strip bcftools plugins/*.so && \
    make DESTDIR=/bcftools_install install && \
    # Building HTSlib explicitly to ensure headers are moved
    make -C htslib-${BCFTOOLS_VERSION} DESTDIR=/bcftools_install install


FROM base

COPY --from=compiler /bcftools_install/usr/local/bin/* /usr/local/bin/
COPY --from=compiler /bcftools_install/usr/local/libexec/bcftools/* /usr/local/libexec/bcftools/

# Install Google Cloud CLI and runtime libraries for BCFtools
RUN apt-get update && \
    # 1. Ensure the keyrings directory exists
    mkdir -p /usr/share/keyrings && \
    # 2. Add Google Cloud repo
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
    > /etc/apt/sources.list.d/google-cloud-sdk.list && \
    # 3. Use curl + gpg with --yes (to overwrite) and --batch
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
    gpg --batch --yes --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    # 4. Final update and install
    apt-get update && \
    apt-get install -y --no-install-recommends \
        google-cloud-cli \
        libcurl4 \
        libgomp1 \
        zlib1g \
        libbz2-1.0 \
        liblzma5 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
