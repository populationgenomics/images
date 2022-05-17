# Add cpg_pipes to driver-base
FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-base:1.1

ENV TZ=Australia \
    DEBIAN_FRONTEND=noninteractive
# hadolint ignore=DL3008, DL3013, DL4006
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    gcc \
    build-essential \
    curl \
    && \
    add-apt-repository ppa:deadsnakes/ppa -y && \
    apt-get install -y --no-install-recommends python3.10 python3.10-distutils python3-apt python3.10-dev  && \
    curl -sS https://bootstrap.pypa.io/get-pip.py | python3.10 && \
    pip3.10 install --no-cache-dir cpg-pipes==0.3.1 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
