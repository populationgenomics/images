FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-r:1.4

# Override this with a specific SHA in images.toml, as the qtl fork doesn't have releases.
ARG VERSION=HEAD

SHELL ["/bin/bash", "-eo", "pipefail", "-c"]

# See https://github.com/saigegit/SAIGE/blob/main/docker/Dockerfile.
# We clone https://github.com/weizhou0/qtl instead of https://github.com/saigegit/SAIGE,
# as the former is an experimental fork for eQTL mapping.
RUN apt-get update && \
    apt-get install -y cmake gfortran libopenblas-base && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    pip3 install cget click && \
    git clone --depth 1 https://github.com/weizhou0/qtl.git && \
    cd qtl && \
    git checkout $VERSION && \
    Rscript extdata/install_packages.R && \
    R CMD INSTALL .

# add function R scripts
RUN mv extdata/step1_fitNULLGLMM_qtl.R extdata/step2_tests_qtl.R \
    extdata/step3_LDmat.R extdata/createSparseGRM.R /usr/local/bin/ && \
    cd .. && \
    rm -r qtl

RUN chmod a+x /usr/local/bin/step1_fitNULLGLMM_qtl.R && \
    chmod a+x /usr/local/bin/step2_tests_qtl.R && \
    chmod a+x /usr/local/bin/step3_LDmat.R && \
    chmod a+x /usr/local/bin/createSparseGRM.R

RUN createSparseGRM.R  --help && \
    step1_fitNULLGLMM_qtl.R --help && \
    step2_tests_qtl.R --help && \
    step3_LDmat.R --help

RUN apt-get update
RUN apt-get install time