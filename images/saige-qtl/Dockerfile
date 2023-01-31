FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-r:1.4

ENV MAMBA_ROOT_PREFIX /root/micromamba

SHELL ["/bin/bash", "-eo", "pipefail", "-c"]
RUN conda create -n RSAIGE r-essentials r-base=4.1.2 python=3.10 \
    conda activate RSAIGE \
    wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir $MAMBA_ROOT_PREFIX && \
    micromamba install -y --prefix $MAMBA_ROOT_PREFIX \
        -c anaconda cmake -c conda-forge gettext lapack r-matrix \
        -c r r-rcpp  r-rcpparmadillo r-data.table r-bh r-matrix \
        -c conda-forge r-spatest r-rcppeigen r-devtools r-skat \
        r-rcppparallel r-optparse boost openblas r-rhpcblasctl r-metaskat \
        r-skat r-qlcmatrix r-rsqlite \
    pip3 install cget click \
    conda env export > environment-RSAIGE.yml