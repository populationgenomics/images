FROM australia-southeast1-docker.pkg.dev/cpg-common/images/nf_bcftools:1.25-1 AS base

ARG VERSION=1.25

# installation code direct from google - https://cloud.google.com/sdk/docs/install#deb
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | \
    tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    wget -O - https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
    gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    apt-get update -y && \
    apt-get install google-cloud-cli --no-install-recommends -y && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*
