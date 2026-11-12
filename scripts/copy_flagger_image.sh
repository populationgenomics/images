#!/usr/bin/env bash
set -ex
# Official HMM-Flagger image (mobinasri/flagger on Docker Hub; no Bioconda package, no quay.io copy) for
# the popgen-pangenome assembly QC stage, run under Nextflow on Seqera. It ships samtools, python3 and
# procps, so it runs as released.
SOURCE_IMAGE="docker://docker.io/mobinasri/flagger:v1.2.0"
IMAGE_NAME="flagger"
IMAGE_TAG="1.2.0"

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy ${SOURCE_IMAGE} docker://australia-southeast1-docker.pkg.dev/cpg-common/images/${IMAGE_NAME}:${IMAGE_TAG}
