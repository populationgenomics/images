FROM python:3.10-bullseye

ARG VERSION=${VERSION:-2.31.0}

RUN apt update && apt install -y \
    apt-transport-https \
    bash \
    ca-certificates \
    curl \
    gnupg \
    wget \
    zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    # Google Cloud SDK: use the script-based installation, as the Debian package is outdated.
    curl https://sdk.cloud.google.com > install.sh && \
    bash install.sh --disable-prompts --install-dir=/opt && \
    rm install.sh && \
    # Download and check the ICA CLI.
    wget -q https://stratus-documentation-us-east-1-public.s3.amazonaws.com/cli/${VERSION}/ica-linux-amd64.sha256 -O /tmp/ica-linux-amd64.sha256 && \
    wget -q https://stratus-documentation-us-east-1-public.s3.amazonaws.com/cli/${VERSION}/ica-linux-amd64.zip -O /tmp/ica-linux-amd64.zip && \
    sed -i 's|/home/ec2-user/workspace/ontrolPlane_ICA_CLI_release_2.31/target/ica-linux-amd64.zip|/tmp/ica-linux-amd64.zip|' /tmp/ica-linux-amd64.sha256 && \
    sha256sum -c /tmp/ica-linux-amd64.sha256 && \
    unzip -d /tmp /tmp/ica-linux-amd64.zip && \
    cp /tmp/linux-amd64/icav2 /usr/local/bin/icav2 && \
    chmod a+x /usr/local/bin/icav2 && \
    rm -rf /tmp/ica-linux-amd64.sha256 /tmp/ica-linux-amd64.zip /tmp/linux-amd64

ENV PATH=$PATH:/opt/google-cloud-sdk/bin
