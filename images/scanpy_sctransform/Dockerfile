
FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-r:1.4

SHELL ["/bin/bash", "-eo", "pipefail", "-c"]

ENV MAMBA_ROOT_PREFIX /root/micromamba
ENV PATH $MAMBA_ROOT_PREFIX/bin:$PATH

RUN apt update && \
    # Some R packages require a C compiler during installation.
    apt install -y build-essential && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    micromamba install -y --prefix $MAMBA_ROOT_PREFIX \
        -c cpg -c bioconda -c conda-forge \
        bioconductor-biomart \
        r-seurat \
        bioconductor-scater  && \
    rm -r /root/micromamba/pkgs && \
    pip3 install --no-cache-dir \
        scanpy==1.9.3 && \
    pip3 install anndata2ri
