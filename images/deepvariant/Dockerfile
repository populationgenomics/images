FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.134.cpg2-1

# Install gsutil if not already present
RUN apt-get update && apt-get install -y \
    curl \
    python3-pip && \
    pip3 install --upgrade pip && \
    pip3 install gsutil

# Set environment variables
ENV DEEPVARIANT_VERSION=1.9.0
ENV DV_DIR=/opt/deepvariant

# Create working directory
WORKDIR ${DV_DIR}

# Copy the binaries from the GCS bucket
RUN gsutil -m cp -r gs://deepvariant/binaries/DeepVariant/${DEEPVARIANT_VERSION}/DeepVariant-${DEEPVARIANT_VERSION}/* ${DV_DIR}

# Make scripts executable
RUN chmod +x ${DV_DIR}/run-prereq.sh ${DV_DIR}/settings.sh

# Run prerequisite and environment setup scripts
RUN ${DV_DIR}/run-prereq.sh && \
    . ${DV_DIR}/settings.sh

# Optionally set entrypoint or CMD
CMD ["/bin/bash"]