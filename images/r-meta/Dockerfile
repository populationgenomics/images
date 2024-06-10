ARG DRIVER_TAG="1.4"

FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-base:$DRIVER_TAG

SHELL ["/bin/bash", "-eo", "pipefail", "-c"]

ENV MAMBA_ROOT_PREFIX /root/micromamba
ENV PATH $MAMBA_ROOT_PREFIX/bin:$PATH

RUN apt update && \
    # Some R packages require a C compiler during installation.
    apt install -y build-essential && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba && \
    mkdir $MAMBA_ROOT_PREFIX && \
    micromamba install -y --prefix $MAMBA_ROOT_PREFIX \
        -c cpg -c bioconda -c conda-forge \
        bioconductor-biomart \
        r-argparser \
        r-arrow \
        r-base \
        r-coloc \
        r-devtools \
        r-essentials \
        r-meta \
        # r-susieR \
        r-tidyverse \
        r-viridis && \
    rm -r /root/micromamba/pkgs && \
    pip3 install --no-cache-dir \
        analysis-runner \
        anndata2ri \
        cloudpathlib[all] \
        cpg-utils \
        cpg-workflows \
        gcsfs \
        pandas
