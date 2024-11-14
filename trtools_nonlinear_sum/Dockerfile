FROM python:3.10-bullseye

ARG VERSION=${VERSION:-piecewise-reg}
ARG GIT_REPO=https://github.com/populationgenomics/TRTools-segreg.git

# Install necessary packages and clean up to reduce image size
RUN apt-get update && apt-get install --no-install-recommends -y \
    git curl bzip2 zip && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/cache/apt/*

# Install Python packages including piecewise_regression
RUN pip install --no-cache-dir metamist piecewise_regression git+${GIT_REPO}@${VERSION}
