#!/usr/bin/env bash
set -ex
# ONT's wf-single-cell image (ontresearch/wf-single-cell on Docker Hub), which runs the single-cell
# steps of popgen-lr-single-cell, ONT's wf-single-cell under Nextflow on Seqera. Multi-tool and built
# by ONT, so it is copied, not rebuilt. Pinned by digest, resolved from the tag
# sha69b55755e85ed2abf209bf19d9f72d7b5dd5f807; skopeo refuses a source with both a tag and a digest,
# so only the digest is given. --all copies every platform, so the copy keeps the source digest. The
# destination tag is the upstream tag, which the pipeline builds from params.wf.container_sha.
SOURCE_IMAGE="docker://docker.io/ontresearch/wf-single-cell@sha256:116d83c0868aa1f366267adcb31ce42a5057d1396fcae7d2069c52ffba6b9cad"
IMAGE_NAME="wf-single-cell"
IMAGE_TAG="sha69b55755e85ed2abf209bf19d9f72d7b5dd5f807"

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy --all ${SOURCE_IMAGE} docker://australia-southeast1-docker.pkg.dev/cpg-common/images/${IMAGE_NAME}:${IMAGE_TAG}
