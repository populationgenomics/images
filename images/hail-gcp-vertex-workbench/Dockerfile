# description of slim base container https://cloud.google.com/vertex-ai/docs/workbench/instances/create-custom-container#slim-base-container
# NB: uses micromamba, NOT conda
FROM us-docker.pkg.dev/deeplearning-platform-release/gcr.io/workbench-container-slim:latest AS base

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    git \
    htop \
    unzip \
    bzip2 \
    zip \
    tar \
    rsync \
    xsltproc pandoc \
    openjdk-11-jdk-headless \
    liblapack3 \
    libopenblas-dev \
    libpq-dev \
    liblz4-dev \
    g++ \
    gcc \
    cmake \
    apt-transport-https \
    ca-certificates \
    gnupg \
    curl && \
    rm -rf /var/lib/apt/lists/*

FROM base AS gcp_build
# gcloud install
RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    apt-get update && \
    apt-get install -y \
    google-cloud-cli

FROM gcp_build AS hail_build
ARG VERSION=${VERSION:-0.2.135}
ENV KERNEL_DISPLAY_NAME="hail-${VERSION}"
ENV ENVIRONMENT_NAME="hail-${VERSION}"
ENV MICROMAMBA_ENV_HOME="/opt/micromamba/envs"
COPY requirements.txt requirements.txt

# create env and install hail + dependences 
RUN micromamba create --prefix "${MICROMAMBA_ENV_HOME}"/"${ENVIRONMENT_NAME}" python=3.10 ipykernel && \
    micromamba run -p "${MICROMAMBA_ENV_HOME}"/"${ENVIRONMENT_NAME}" pip install --no-cache-dir --no-deps -r requirements.txt && \
    micromamba run -p "${MICROMAMBA_ENV_HOME}"/"${ENVIRONMENT_NAME}" python -m ipykernel install --name "${ENVIRONMENT_NAME}" --display-name "${KERNEL_DISPLAY_NAME}" && \
    curl -o "$(micromamba run -p "${MICROMAMBA_ENV_HOME}"/"${ENVIRONMENT_NAME}" find_spark_home.py)/jars/gcs-connector-hadoop2-2.2.28.jar" https://storage.googleapis.com/hadoop-lib/gcs/gcs-connector-hadoop2-2.2.28.jar


