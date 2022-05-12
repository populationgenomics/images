# Standard analysis-runner driver image with cpg_pipes. 
FROM ubuntu:20.04

RUN apt update && apt install -y \
        apt-transport-https \
        bash \
        bzip2 \
        ca-certificates \
        curl \
        git \
        gnupg \
        jq \
        wget \
        python3 \
        python3-pip \
        zip && \
    pip install --no-cache-dir cpg_pipes  && \
    # Google Cloud SDK
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt update && apt install -y google-cloud-sdk && \
    # skopeo
    . /etc/os-release && \
    echo "deb https://download.opensuse.org/repositories/devel:/kubic:/libcontainers:/stable/xUbuntu_${VERSION_ID}/ /" | tee /etc/apt/sources.list.d/devel:kubic:libcontainers:stable.list && \
    curl https://download.opensuse.org/repositories/devel:/kubic:/libcontainers:/stable/xUbuntu_${VERSION_ID}/Release.key | apt-key add - && \
    apt-get update && \
    apt-get -y install skopeo && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*