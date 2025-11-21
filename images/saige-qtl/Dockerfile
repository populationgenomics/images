FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-r:1.4

# Explicitly specify desired image label and corresponding commit, as the qtl fork doesn't have releases.
# Currently building from the main_initVal branch
ARG VERSION=0.3.2.1-c2c9fb2
ARG VERSION_COMMIT=c2c9fb20a0efc4c9a4e7289f6f444d3e19f32ec9

SHELL ["/bin/bash", "-eo", "pipefail", "-c"]

# See https://github.com/saigegit/SAIGE/blob/main/docker/Dockerfile.
# We clone https://github.com/weizhou0/qtl instead of https://github.com/saigegit/SAIGE,
# as the former is an experimental fork for eQTL mapping.
RUN apt-get update && \
    apt-get install -y cmake gfortran libopenblas-base && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    micromamba install -y --prefix $MAMBA_ROOT_PREFIX -c bioconda -c conda-forge libgcc-ng==12.2.0 savvy superlu && \
    rm -r $MAMBA_ROOT_PREFIX/pkgs && \
    pip3 install cget click && \
    git clone https://github.com/weizhou0/qtl.git && \
    cd qtl && \
    git checkout $VERSION_COMMIT && \
    sed -i $'/#include.*cassert/a\\\n#include <tbb/concurrent_vector.h>\n' src/GENO_null.hpp && \
    Rscript extdata/install_packages.R && \
    R CMD INSTALL . && \
    for f in step1_fitNULLGLMM_qtl.R step2_tests_qtl.R step3_gene_pvalue_qtl.R makeGroupFile.R; do \
        sed '/#!.*env.*pixi.*Rscript$/s,/.*,'$MAMBA_ROOT_PREFIX'/bin/Rscript,' < extdata/$f > /usr/local/bin/$f; \
        chmod a+x /usr/local/bin/$f; \
    done && \
    mv extdata/input/n.indep_100_n.cell_1_01.step1* \
       extdata/input/seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt \
       extdata/input/genotype_10markers.vcf.gz \
       extdata/input/genotype_10markers.vcf.gz.csi \
       extdata/input/gene_1_cis_region.txt /usr/local/bin/ && \
    cd .. && \
    rm -r qtl
