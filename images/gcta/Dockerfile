FROM python:3.10-bullseye

ARG VERSION=${VERSION:-1.94.1}

ENV MAMBA_ROOT_PREFIX /root/micromamba
ENV PATH $MAMBA_ROOT_PREFIX/bin:$PATH

RUN apt-get update && apt-get install -y git wget bash bzip2 zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir ${MAMBA_ROOT_PREFIX} && \
    micromamba install -y --prefix ${MAMBA_ROOT_PREFIX} -c bioconda -c conda-forge \
    gcta=${VERSION} \
    google-cloud-sdk \
    && \
    rm -r /root/micromamba/pkgs
