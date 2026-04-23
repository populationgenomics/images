FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.137.cpg1-2 AS base

ARG VERSION=${VERSION:-2.0.3}

RUN git clone https://github.com/hassansaei/vntyper.git && \
    cd vntyper && \
    git checkout v${VERSION} && \
    pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir . && \
    vntyper install-references -d /vntyper/references --references hg38
