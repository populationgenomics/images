FROM python:3.10-bullseye

ARG VERSION=${VERSION:-largest_allele_effect}
ARG GIT_REPO=https://github.com/populationgenomics/TRTools-nonlinear-sum.git

# Install necessary packages and clean up to reduce image size
RUN apt-get update && apt-get install --no-install-recommends -y \
    git curl bzip2 zip && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/cache/apt/*

RUN pip install --no-cache-dir metamist git+${GIT_REPO}@${VERSION}