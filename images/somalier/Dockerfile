FROM python:3.11-slim-trixie

ENV DEBIAN_FRONTEND=noninteractive

ARG VERSION=v0.3.6

RUN apt-get update && apt-get install -y \
    bc \
    bzip2 \
    curl \
    git \
    gnupg \
    procps \
    wget \
    zip  \
    && \
    rm -rf /var/lib/apt/lists/* &&  \
    rm -rf /tmp/*

# Install Google Cloud SDK (gcloud)
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | \
    tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    wget -O - https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
    gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    apt-get update -y && \
    apt-get install google-cloud-cli --no-install-recommends -y && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*
# Install somalier
RUN wget https://github.com/brentp/somalier/releases/download/v${VERSION}/somalier && \
    mv ./somalier /usr/local/bin/somalier && \
    chmod +x /usr/local/bin/somalier
