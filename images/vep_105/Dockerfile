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
    echo $VEP_SHARE && \
    rm $VEP_SHARE/TissueExpression.pm \
        $VEP_SHARE/ancestral.pm \
        $VEP_SHARE/context.pm \
        $VEP_SHARE/de_novo_donor.pl \
        $VEP_SHARE/extended_splice.pl \
        $VEP_SHARE/gerp_dist.pl \
        $VEP_SHARE/loftee_splice_utils.pl \
        $VEP_SHARE/splice_site_scan.pl \
        $VEP_SHARE/svm.pl \
        $VEP_SHARE/utr_splice.pl && \
    git clone https://github.com/populationgenomics/loftee_38.git && \
    cp -r loftee_38/* $VEP_SHARE && \
    rm -rf loftee_38 && \
    # Install gcsfuse to mount VEP cache files
    GCSFUSE_REPO=gcsfuse-$(lsb_release -c -s) && \
    echo "deb https://packages.cloud.google.com/apt $GCSFUSE_REPO main" | tee /etc/apt/sources.list.d/gcsfuse.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && \
    apt-get install -y gcsfuse && \
    # Install Bio::Perl for Loftee (conda's perl-bioperl doesn't work)
    # run this last, as other commands seem to re-break the installation
    vep_install --AUTO a --NO_UPDATE --NO_HTSLIB && \
    ln -fs $MAMBA_ROOT_PREFIX/share/ensembl-vep-$VERSION* $MAMBA_ROOT_PREFIX/share/ensembl-vep && \
    echo "export PERL5LIB=$MAMBA_ROOT_PREFIX/share/ensembl-vep" >> ~/.bashrc && \
    echo "export LOFTEE_PLUGIN_PATH=$MAMBA_ROOT_PREFIX/share/ensembl-vep" >> ~/.bashrc
