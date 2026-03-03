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
    # UNSET conflicting paths so the compiler looks in /usr/include
    unset CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH LIBRARY_PATH && \
    ./configure --enable-libcurl --enable-s3 --enable-gcs && \
    make && \
    strip bcftools plugins/*.so && \
    make DESTDIR=/bcftools_install install && \
    # The internal htslib directory name in the release tarball 
    # is usually just 'htslib' or 'htslib-<version>'
    make -C htslib-${BCFTOOLS_VERSION} DESTDIR=/bcftools_install install


FROM base

COPY --from=compiler /bcftools_install/usr/local/bin/* /usr/local/bin/
COPY --from=compiler /bcftools_install/usr/local/libexec/bcftools/* /usr/local/libexec/bcftools/

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | \
    tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    wget -O - https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
    gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    apt-get update -y && \
    apt-get install google-cloud-cli --no-install-recommends -y && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*
