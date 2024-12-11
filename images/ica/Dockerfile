FROM debian:buster-slim

ARG VERSION=${VERSION:-2.31.0}

RUN apt-get update && apt-get install -y git wget bash bzip2 zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    wget -q https://stratus-documentation-us-east-1-public.s3.amazonaws.com/cli/${VERSION}/ica-linux-amd64.sha256 -O /tmp/ica-linux-amd64.sha256 && \
    wget -q https://stratus-documentation-us-east-1-public.s3.amazonaws.com/cli/${VERSION}/ica-linux-amd64.zip -O /tmp/ica-linux-amd64.zip && \
    sed -i 's|/home/ec2-user/workspace/ontrolPlane_ICA_CLI_release_2.31/target/ica-linux-amd64.zip|/tmp/ica-linux-amd64.zip|' /tmp/ica-linux-amd64.sha256 \
    sha256sum -c /tmp/ica-linux-amd64.sha256 && \
    unzip -d /tmp /tmp/ica-linux-amd64.zip && \
    cp /tmp/linux-amd64/icav2 /usr/local/bin/icav2 && \
    chmod a+x /usr/local/bin/icav2 && \
    rm -rf /tmp/ica-linux-amd64.sha256 /tmp/ica-linux-amd64.zip /tmp/linux-amd64