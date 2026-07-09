# Use a slim Python image for efficiency
FROM python:3.11-slim

ENV VERSION=0.0.1

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && \
    apt-get install -y \
      git \
      procps && \
    rm -rf /var/lib/apt/lists/*

# (Nextflow's bin/ script is automatically added to the PATH by Nextflow,
#  so we don't strictly need to COPY it into the container if we mount it,
#  but for a truly self-contained container on GCP AR, we can include it.)
RUN pip install --no-cache-dir \
    "google-cloud-dataproc>=5.0.0" \
    "google-cloud-storage>=2.0.0" \
    "metamist"

# install the latest cpg-utils code for testing the dataproc wrapper and analysis registration \
RUN git clone -b add_dataproc_handler https://github.com/populationgenomics/cpg-utils.git && \
    cd cpg-utils && \
    pip install .
