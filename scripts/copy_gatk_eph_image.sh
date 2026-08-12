#!/usr/bin/env bash
set -ex
SOURCE_IMAGE="docker://us.gcr.io/broad-dsde-methods/eph/gatk:eph_federate_svs-53b4f16a7"
IMAGE_NAME="gatk"
IMAGE_TAG="eph_federate_svs-53b4f16a7"

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy ${SOURCE_IMAGE} docker://australia-southeast1-docker.pkg.dev/cpg-common/images/${IMAGE_NAME}:${IMAGE_TAG}
