FROM python:3.12-bullseye

ARG VERSION=${VERSION:-6.1.2}

RUN apt update && apt install -y --no-install-recommends \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        cmake \
        gfortran \
        git \
        gnupg \
        libopenblas-dev \
        openjdk-11-jdk-headless \
        wget \
        zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

# install a skinny version of gcloud
# gcloud is a monster - contains thousands of man pages, a whole python installation, etc.0
# gcloud-lite is a much smaller version, with only the bits we need
# gcloud-go would be even smaller, but I couldn't grok the install
ARG CLOUD_SDK_VERSION=499.0.0
ENV CLOUD_SDK_VERSION=$CLOUD_SDK_VERSION
ENV CLOUDSDK_PYTHON="/usr/local/bin/python"
ADD "https://github.com/tonymet/gcloud-lite/releases/download/${CLOUD_SDK_VERSION}/google-cloud-cli-${CLOUD_SDK_VERSION}-linux-x86_64-lite.tar.gz" gcloud.tar.gz
RUN tar -zxf gcloud.tar.gz \
    && mv google-cloud-sdk /bin/google-cloud-sdk \
    && rm gcloud.tar.gz \
    && gcloud components reinstall

ENV PATH="/bin/google-cloud-sdk/bin:$PATH"

# install nextflow
ADD https://get.nextflow.io nextflow
RUN chmod +x nextflow && \
    mv nextflow /usr/bin && \
    nextflow self-update

# botocore is required for the Gcloud CLI
RUN pip install --no-cache-dir metamist git+https://github.com/populationgenomics/talos.git@${VERSION}
