#!/usr/bin/env bash
set -ex
# Google's DeepVariant, released 1.10.0 (tag 1.10.0), pinned by digest. popgen-pangenome's assembly QC
# runs it on HiFi reads aligned to each sample's own assembly. skopeo refuses a source with both a tag
# and a digest, so only the digest is given. The tag is a single linux/amd64 image, not a platform list,
# so the copy keeps this digest. The pangenome-aware release sits beside it as
# deepvariant:pangenome_aware_deepvariant-1.10.0 (copy_deepvariant_pangenome_image.sh).
SOURCE_IMAGE="docker://docker.io/google/deepvariant@sha256:962e5a83b1d76aae6990625d47102785f791603f2138aa1fa9aa4fb6a2eecbe6"
IMAGE_NAME="deepvariant"
IMAGE_TAG="1.10.0"

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy ${SOURCE_IMAGE} docker://australia-southeast1-docker.pkg.dev/cpg-common/images/${IMAGE_NAME}:${IMAGE_TAG}
