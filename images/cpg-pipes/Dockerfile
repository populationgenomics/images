# Add cpg_pipes to driver-base
FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-base:1.1
# hadolint ignore=DL3008, DL3013
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    build-essential \
    python3 \
    python3-pip \
    python3-dev \
    && \
    pip install --no-cache-dir cpg-pipes==0.3.1 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
