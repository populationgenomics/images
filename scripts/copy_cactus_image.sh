#!/usr/bin/env bash
set -ex
# Official cactus image for running the pangenome pipeline
SOURCE_IMAGE="docker://quay.io/comparative-genomics-toolkit/cactus:v3.3.0"
IMAGE_NAME="cactus"
IMAGE_TAG="3.3.0"

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy ${SOURCE_IMAGE} docker://australia-southeast1-docker.pkg.dev/cpg-common/images/${IMAGE_NAME}:${IMAGE_TAG}
