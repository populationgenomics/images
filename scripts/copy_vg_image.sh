#!/usr/bin/env bash
set -ex
# Official vg 1.77.0 image (tag v1.77.0), pinned by digest. skopeo refuses a source with both a tag
# and a digest, so only the digest is given. The tag is a multi-platform list, so --all copies every
# platform and the copy keeps the source digest.
SOURCE_IMAGE="docker://quay.io/vgteam/vg@sha256:31edee2fa9386d5c47088d303a97a14a803f7446a6d9b0bdc469310a1b94da5c"
IMAGE_NAME="vg"
IMAGE_TAG="1.77.0"

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy --all ${SOURCE_IMAGE} docker://australia-southeast1-docker.pkg.dev/cpg-common/images/${IMAGE_NAME}:${IMAGE_TAG}
