FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.134.cpg2-1

ARG DV_GPU_BUILD=0
ARG VERSION=1.9.0
ARG PYTHON_VERSION=3.10
ARG TF_ENABLE_ONEDNN_OPTS=1

ENV DV_GPU_BUILD=${DV_GPU_BUILD}
ENV VERSION=${VERSION}
ENV PYTHON_VERSION=${PYTHON_VERSION}
ENV TF_ENABLE_ONEDNN_OPTS=${TF_ENABLE_ONEDNN_OPTS}
ENV DV_DIR=/opt/deepvariant
ENV DV_BIN_PATH=/opt/deepvariant/bin

WORKDIR ${DV_DIR}

RUN gsutil -m cp -r gs://deepvariant/binaries/DeepVariant/${VERSION}/DeepVariant-${VERSION}/* ${DV_DIR}
RUN gsutil cp gs://cpg-bioheart-test/deepvariant_test/run_deepvariant ${DV_BIN_PATH}/run_deepvariant && \
    gsutil cp gs://cpg-bioheart-test/deepvariant_test/run_deepvariant.py ${DV_BIN_PATH}/run_deepvariant.py && \
    chmod +x ${DV_BIN_PATH}/run_deepvariant

RUN chmod +x ${DV_DIR}/run-prereq.sh ${DV_DIR}/settings.sh

RUN ${DV_DIR}/run-prereq.sh && \
    bash -c "source ${DV_DIR}/settings.sh"

CMD ["/bin/bash"]
