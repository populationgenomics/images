FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.137.cpg1-2 AS base

ARG VERSION=${VERSION:-2.0.3}

RUN conda install -y -c bioconda -c conda-forge \
        bwa=0.7.18 \
        samtools=1.20 \
        fastp=0.23.4 \
    && conda clean -afy

RUN git clone https://github.com/hassansaei/vntyper.git && \
    cd vntyper && \
    git checkout v${VERSION} && \
    pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir . && \
    mkdir -p /vntyper/references && \
    vntyper install-references --output-dir /vntyper/references --references hg38
