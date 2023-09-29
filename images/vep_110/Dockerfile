FROM debian:bullseye-slim

ENV MAMBA_ROOT_PREFIX=/root/micromamba
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH
ARG VERSION=${VERSION:-110.1}

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
        zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir $MAMBA_ROOT_PREFIX && \
    micromamba install -y --prefix $MAMBA_ROOT_PREFIX -c bioconda -c conda-forge \
        ensembl-vep=${VERSION} \
        google-cloud-sdk \
        # Loftee deps:
        perl-bio-bigfile \
        perl-dbd-sqlite \
        perl-list-moreutils \
        samtools && \
    rm -r /root/micromamba/pkgs && \
    # remove broken (GRCh37) LOFTEE files, copy in new GRCh38 content
    VEP_SHARE=$MAMBA_ROOT_PREFIX/share/ensembl-vep-$VERSION-* && \
    git clone -b grch38 https://github.com/konradjk/loftee.git && \
    mv loftee/* $VEP_SHARE && \
    rm -rf loftee && \
    # pull AlphaMissense \
    curl https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/110/AlphaMissense.pm > ${VEP_SHARE}/AlphaMissense.pm \
    # Install gcsfuse to mount VEP cache files
    GCSFUSE_REPO=gcsfuse-$(lsb_release -c -s) && \
    echo "deb http://packages.cloud.google.com/apt $GCSFUSE_REPO main" | tee /etc/apt/sources.list.d/gcsfuse.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && \
    apt-get install -y gcsfuse && \
    # Install Bio::Perl for Loftee (conda's perl-bioperl doesn't work)
    # run this last, as other commands seem to re-break the installation
    vep_install --AUTO a --NO_UPDATE --NO_HTSLIB && \
    ln -fs $VEP_SHARE $MAMBA_ROOT_PREFIX/share/ensembl-vep && \
    echo "export PERL5LIB=$MAMBA_ROOT_PREFIX/share/ensembl-vep" >> ~/.bashrc && \
    echo "export LOFTEE_PLUGIN_PATH=$MAMBA_ROOT_PREFIX/share/ensembl-vep" >> ~/.bashrc
