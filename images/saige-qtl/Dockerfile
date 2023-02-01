FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-r:1.4

SHELL ["/bin/bash", "-eo", "pipefail", "-c"]
RUN micromamba install -y --prefix $MAMBA_ROOT_PREFIX \
        -c cpg -c bioconda -c conda-forge -c anaconda -c r \
        cmake gettext lapack r-matrix \
        r-rcpp r-rcpparmadillo r-data.table r-bh r-matrix \
        r-spatest r-rcppeigen r-devtools r-skat \
        r-rcppparallel r-optparse boost openblas r-rhpcblasctl r-metaskat \
        r-skat r-qlcmatrix r-rsqlite && \
    rm -r /root/micromamba/pkgs && \
    pip3 install cget click

RUN FLAGPATH=`which python | sed 's|/bin/python$||'` && \
    export LDFLAGS="-L${FLAGPATH}/lib" && \
    export CPPFLAGS="-I${FLAGPATH}/include"

RUN src_branch=main && \
    repo_src_url=https://github.com/weizhou0/qtl && \
    git clone --depth 1 -b $src_branch $repo_src_url && \
    Rscript ./qtl/extdata/install_packages.R && \
    R CMD INSTALL --library=/Users/anncuo/Documents/Garvan/software/SAIGE SAIGE