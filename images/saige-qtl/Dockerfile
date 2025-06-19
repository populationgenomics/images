FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-r:1.4

# Override this with a specific SHA in images.toml, as the qtl fork doesn't have releases.
ARG VERSION=c3ddadda416a0be349872ed3b486ec020c371764

SHELL ["/bin/bash", "-eo", "pipefail", "-c"]

# See https://github.com/saigegit/SAIGE/blob/main/docker/Dockerfile.
# We clone https://github.com/weizhou0/qtl instead of https://github.com/saigegit/SAIGE,
# as the former is an experimental fork for eQTL mapping.
RUN apt-get update && \
    apt-get install -y cmake gfortran libopenblas-base && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    pip3 install cget click && \
    git clone https://github.com/weizhou0/qtl.git && \
    cd qtl && \
    git checkout $VERSION && \
    Rscript extdata/install_packages.R && \
    R CMD INSTALL . && \
    mv extdata/step1_fitNULLGLMM_qtl.R \
       extdata/step2_tests_qtl.R \
       extdata/step3_gene_pvalue_qtl.R \
       extdata/input/n.indep_100_n.cell_1_01.step1* \
       extdata/input/seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt \
       extdata/input/genotype_10markers.vcf.gz \
       extdata/input/genotype_10markers.vcf.gz.csi \
       extdata/input/gene_1_cis_region.txt /usr/local/bin/ && \
    chmod a+x /usr/local/bin/step1_fitNULLGLMM_qtl.R \
              /usr/local/bin/step2_tests_qtl.R \
              /usr/local/bin/step3_gene_pvalue_qtl.R && \
    cd .. && \
    rm -r qtl
