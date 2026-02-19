FROM debian:trixie-slim AS base

ENV DEBIAN_FRONTEND=noninteractive
ARG VERSION=${VERSION:-0.3.1}

# Install dependencies
RUN apt-get update && apt-get install -y \
    git wget bzip2 zip curl bc && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*
# Install Google Cloud SDK (gcloud)
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | \
    tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
    apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update && apt-get install -y google-cloud-cli && \
    rm -r /var/lib/apt/lists/*
# Install somalier
RUN wget https://github.com/brentp/somalier/releases/download/v${VERSION}/somalier && \
    mv ./somalier /usr/local/bin/somalier && \
    chmod +x /usr/local/bin/somalier
