#!/usr/bin/env bash
set -ex
# Google's pangenome-aware DeepVariant, released 1.10.0 (tag pangenome_aware_deepvariant-1.10.0),
# pinned by digest. skopeo refuses a source with both a tag and a digest, so only the digest is given.
SOURCE_IMAGE="docker://docker.io/google/deepvariant@sha256:996e456cd839e29abecd0fd3a9f552e8e09b59f945acfeb13022a2b88f9a2ede"
IMAGE_NAME="deepvariant"
IMAGE_TAG="pangenome_aware_deepvariant-1.10.0"

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy ${SOURCE_IMAGE} docker://australia-southeast1-docker.pkg.dev/cpg-common/images/${IMAGE_NAME}:${IMAGE_TAG}
