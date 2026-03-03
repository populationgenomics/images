#!/usr/bin/env bash
set -ex
SOURCE_IMAGE="docker://broadinstitute/gatk-nightly:2026-03-01-4.6.2.0-21-ge8c49f600-NIGHTLY-SNAPSHOT"
IMAGE_NAME="gatk"
IMAGE_TAG="4.6.2.0-21-ge8c49f600-NIGHTLY-SNAPSHOT"

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy ${SOURCE_IMAGE} docker://australia-southeast1-docker.pkg.dev/cpg-common/images/${IMAGE_NAME}:${IMAGE_TAG}
