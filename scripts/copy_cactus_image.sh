#!/usr/bin/env bash
set -ex
# Official cactus image for running tha pangenome pipeline
SOURCE_IMAGE="quay.io/comparative-genomics-toolkit/cactus:3.3.0"
IMAGE_NAME="cactus"
IMAGE_TAG="3.3.0"

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy ${SOURCE_IMAGE} docker://australia-southeast1-docker.pkg.dev/cpg-common/images/${IMAGE_NAME}:${IMAGE_TAG}
