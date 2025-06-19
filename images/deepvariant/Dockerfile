FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.134.cpg2-1

# Install gsutil if not already present
RUN apt-get update && apt-get install -y \
    curl \
    python3-pip && \
    pip3 install --upgrade pip && \
    pip3 install gsutil

ENV DEEPVARIANT_VERSION=1.9.0
ENV DV_DIR=/opt/deepvariant

WORKDIR ${DV_DIR}

RUN gsutil -m cp -r gs://deepvariant/binaries/DeepVariant/${DEEPVARIANT_VERSION}/DeepVariant-${DEEPVARIANT_VERSION}/* ${DV_DIR}

RUN chmod +x ${DV_DIR}/run-prereq.sh ${DV_DIR}/settings.sh

RUN ${DV_DIR}/run-prereq.sh && \
    . ${DV_DIR}/settings.sh

CMD ["/bin/bash"]