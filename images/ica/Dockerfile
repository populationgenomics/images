FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:latest

ARG VERSION=${VERSION:-2.31.0}

# Download and check the ICA CLI.
RUN wget -q https://stratus-documentation-us-east-1-public.s3.amazonaws.com/cli/${VERSION}/ica-linux-amd64.sha256 -O /tmp/ica-linux-amd64.sha256 && \
    wget -q https://stratus-documentation-us-east-1-public.s3.amazonaws.com/cli/${VERSION}/ica-linux-amd64.zip -O /tmp/ica-linux-amd64.zip && \
    sed -i 's|/home/ec2-user/workspace/ontrolPlane_ICA_CLI_release_2.31/target/ica-linux-amd64.zip|/tmp/ica-linux-amd64.zip|' /tmp/ica-linux-amd64.sha256 && \
    sha256sum -c /tmp/ica-linux-amd64.sha256 && \
    unzip -d /tmp /tmp/ica-linux-amd64.zip && \
    cp /tmp/linux-amd64/icav2 /usr/local/bin/icav2 && \
    chmod a+x /usr/local/bin/icav2 && \
    rm -rf /tmp/ica-linux-amd64.sha256 /tmp/ica-linux-amd64.zip /tmp/linux-amd64

RUN pip install git+https://github.com/Illumina/ica-sdk-python.git \
    && pip install typing-extensions --upgrade

ENV PATH=$PATH:/opt/google-cloud-sdk/bin
