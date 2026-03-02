ARG VERSION=4.6.2.0

FROM python:3.10-slim-bookworm AS base

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        curl \
        git \
        git-lfs \
        gnupg \
        jq \
        libgomp1 \
        openjdk-17-jdk-headless \
        procps \
        rsync \
        software-properties-common \
        unzip \
        wget \
        zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    pip install --no-cache-dir --upgrade pip

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
    ./configure --enable-libcurl --enable-s3 --enable-gcs && \
    make && \
    strip bcftools plugins/*.so && \
    make DESTDIR=/bcftools_install install && \
    make -C htslib-$BCFTOOLS_VERSION DESTDIR=/bcftools_install install


# build GATK image from source, using the compiled bcftools binaries and plugins
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

# Clone GATK source code (including large files with git-lfs) and build the GATK jar file
RUN git clone --branch ${VERSION} --depth 1 https://github.com/broadinstitute/gatk.git

RUN git lfs install && \
    git lfs pull
# Build GATK with htsjdk 4.3.0 for CRAM v3.1 support
RUN ./gradlew clean printVersion shadowJar -Dhtsjdk.version=4.3.0

ENV PATH=$PATH:/gatk-4.6.2.0
WORKDIR /gatk-4.6.2.0

# Install Miniconda3-py310_23.10.0-1
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.10.0-1-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -p /opt/miniconda -b && \
    rm /tmp/miniconda.sh && \
    /opt/miniconda/bin/conda clean -afy && \
    conda config --set auto_update_conda false && \
    conda config --set solver libmamba && \
    ./gradlew localDevCondaEnv && \
    source activate gatk
