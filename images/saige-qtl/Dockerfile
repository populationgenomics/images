FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-r:1.4

SHELL ["/bin/bash", "-eo", "pipefail", "-c"]
RUN micromamba install -y --prefix $MAMBA_ROOT_PREFIX \
        -c cpg -c bioconda -c conda-forge -c anaconda \
        r-essentials r-base=4.1.2 python=3.10 \
        cmake gettext lapack r-matrix \
        r r-rcpp r-rcpparmadillo r-data.table r-bh r-matrix \
        r-spatest r-rcppeigen r-devtools r-skat \
        r-rcppparallel r-optparse boost openblas r-rhpcblasctl r-metaskat \
        r-skat r-qlcmatrix r-rsqlite && \
    rm -r /root/micromamba/pkgs && \
    pip3 install cget click