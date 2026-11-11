#!/usr/bin/env bash
set -ex
# ONT's shared tools image (ontresearch/wf-common on Docker Hub), which runs ingress (fastcat), BAM
# merging and the report in popgen-lr-single-cell, ONT's wf-single-cell under Nextflow on Seqera.
# Multi-tool and built by ONT, so it is copied, not rebuilt. Pinned by digest, resolved from the tag
# shafdd79f8e4a6faad77513c36f623693977b92b08e; skopeo refuses a source with both a tag and a digest,
# so only the digest is given. The tag is a multi-platform list (amd64, arm64), so --all copies every
# platform and the copy keeps the source digest. The destination tag is the upstream tag, which the
# pipeline builds from params.wf.common_sha.
SOURCE_IMAGE="docker://docker.io/ontresearch/wf-common@sha256:fd1faa356004d0f4062625fef50782cef0865fd0c78e050c7f875a39f32410a8"
IMAGE_NAME="wf-common"
IMAGE_TAG="shafdd79f8e4a6faad77513c36f623693977b92b08e"

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy --all ${SOURCE_IMAGE} docker://australia-southeast1-docker.pkg.dev/cpg-common/images/${IMAGE_NAME}:${IMAGE_TAG}
