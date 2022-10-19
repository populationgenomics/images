FROM debian:stretch-slim

ENV MAMBA_ROOT_PREFIX=/root/micromamba
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH
ARG VERSION=${VERSION:-105.0}

# hadolint ignore=DL3008
RUN apt-get update && \
    apt-get install -y \
    bash \
    bzip2 \
    curl \
    git \
    gnupg2 \
    lsb-release \
    wget \
    zip \
    && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir $MAMBA_ROOT_PREFIX && \
    micromamba install -y --prefix $MAMBA_ROOT_PREFIX -c bioconda -c conda-forge \
    ensembl-vep=${VERSION} \
    # Loftee deps:
    perl-dbd-sqlite \
    perl-list-moreutils \
    samtools \
    google-cloud-sdk \
    && \
    rm -r /root/micromamba/pkgs

# Install Bio::Perl for Loftee (conda's perl-bioperl doesn't work)
RUN vep_install --AUTO a --NO_HTSLIB

# Install gcsfuse to mount VEP cache files
RUN GCSFUSE_REPO=gcsfuse-$(lsb_release -c -s) && \
    echo "deb http://packages.cloud.google.com/apt $GCSFUSE_REPO main" | tee /etc/apt/sources.list.d/gcsfuse.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install --no-install-recommends -y gcsfuse

RUN git clone -b grch38 https://github.com/konradjk/loftee.git && \
    cp -r loftee/* $MAMBA_ROOT_PREFIX/share/ensembl-vep-${VERSION}* \
    && rm -rf loftee
