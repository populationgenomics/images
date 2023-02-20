# Similar to analysis-runner.docker.base without gcloud SDK
# https://github.com/populationgenomics/analysis-runner/blob/main/driver/Dockerfile.base
# https://github.com/populationgenomics/analysis-runner/blob/main/driver/Dockerfile.hail
FROM python:3.10-bullseye

RUN apt update && apt install -y \
        apt-transport-https \
        bash \
        build-essential \
        bzip2 \
        ca-certificates \
        curl \
        g++ \
        gcc \
        git \
        gnupg \
        jq \
        libfontconfig \
        liblapack3 \
        libopenblas-base \
        openjdk-11-jdk-headless \
        python3 \
        rsync \
        software-properties-common \
        wget \
        zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* \
    cd /opt && \
    curl -sSL https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-linux-x86_64.tar.bz2 | tar xj && \
    cd .. && \
    printf '#!/bin/bash\nexport OPENSSL_CONF=/dev/null\n/opt/phantomjs-2.1.1-linux-x86_64/bin/phantomjs "$@"\n' > /usr/local/bin/phantomjs && \
    chmod +x /usr/local/bin/phantomjs && \
    # Install Hail from the CPG fork.
    git clone https://github.com/populationgenomics/hail.git && \
    cd hail && \
    git checkout $HAIL_SHA && \
    cd hail && \
    # Install locally, avoiding the need for a pip package.
    # DEPLOY_REMOTE avoids a dev suffix being appended to dataproc initialization paths.
    make install DEPLOY_REMOTE=1 && \
    cd ../.. && \
    rm -rf hail && \
    pip --no-cache-dir install \
        cloudpathlib[all] \
        cpg-utils \
        cpg-workflows \
        gcsfs \
        hail \
        sample-metadata \
        statsmodels

