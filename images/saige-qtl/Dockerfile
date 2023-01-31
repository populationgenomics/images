FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-r:1.4

ENV MAMBA_ROOT_PREFIX /root/micromamba

SHELL ["/bin/bash", "-eo", "pipefail", "-c"]
RUN wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir $MAMBA_ROOT_PREFIX && \
    micromamba install -y --prefix $MAMBA_ROOT_PREFIX \
        -c cpg -c bioconda -c conda-forge \
        r-essentials r-base=4.1.2 python=3.10 \
        anaconda cmake \
        -c conda-forge gettext lapack r-matrix \
        r r-rcpp r-rcpparmadillo r-data.table r-bh r-matrix \
        r-spatest r-rcppeigen r-devtools r-skat \
        r-rcppparallel r-optparse boost openblas r-rhpcblasctl r-metaskat \
        r-skat r-qlcmatrix r-rsqlite \
    rm -r /root/micromamba/pkgs \
    pip3 install cget click