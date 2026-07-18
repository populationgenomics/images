FROM debian:buster-slim

ENV MAMBA_ROOT_PREFIX /root/micromamba
ENV PATH $MAMBA_ROOT_PREFIX/bin:$PATH
ARG VERSION=${VERSION:-2.28}

RUN apt-get update && apt-get install -y git wget bash bzip2 zip curl bc && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir ${MAMBA_ROOT_PREFIX} && \
    micromamba install -y --prefix ${MAMBA_ROOT_PREFIX} -c bioconda -c conda-forge \
    samtools google-cloud-sdk && \
    rm -r /root/micromamba/pkgs && \
    curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -jxvf - && \
    mv ./minimap2-2.28_x64-linux/minimap2 /root/micromamba/bin/minimap2 && \
    chmod +x /root/micromamba/bin/minimap2 && \
    rm -rf ./minimap2-2.28_x64-linux
