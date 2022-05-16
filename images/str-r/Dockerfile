FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-r:1.1

ENV MAMBA_ROOT_PREFIX /root/micromamba
ENV PATH $MAMBA_ROOT_PREFIX/bin:$PATH

SHELL ["/bin/bash", "-eo", "pipefail", "-c"]
RUN wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir $MAMBA_ROOT_PREFIX && \
    micromamba install -y --prefix $MAMBA_ROOT_PREFIX \
        -c cpg -c bioconda -c conda-forge \
        r-argparser \
        r-base=4.1.1 \
        r-essentials \
        r-tidyverse \
        r-googleCloudStorageR \
        r-hrbrthemes \
        r-vcfR \
        r-readr \
        r-kableExtra \
        r-gridExtra \
        r-DT \
        r-plotly && \
    rm -r /root/micromamba/pkgs
