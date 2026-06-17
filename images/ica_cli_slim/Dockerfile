FROM google/cloud-sdk:slim

ARG VERSION=${VERSION:-2.31.0}

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget unzip ca-certificates && \
    wget -q https://stratus-documentation-us-east-1-public.s3.amazonaws.com/cli/${VERSION}/ica-linux-amd64.sha256 -O /tmp/ica.sha256 && \
    wget -q https://stratus-documentation-us-east-1-public.s3.amazonaws.com/cli/${VERSION}/ica-linux-amd64.zip -O /tmp/ica.zip && \
    echo "$(cut -d' ' -f1 /tmp/ica.sha256)  /tmp/ica.zip" | sha256sum -c - && \
    unzip -d /tmp /tmp/ica.zip && \
    install -m 0755 /tmp/linux-amd64/icav2 /usr/local/bin/icav2 && \
    apt-get purge -y wget unzip && apt-get autoremove -y && \
    rm -rf /tmp/ica.* /tmp/linux-amd64 /var/lib/apt/lists/*
