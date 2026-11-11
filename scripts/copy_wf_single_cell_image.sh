#!/usr/bin/env bash
set -ex
# ONT's wf-single-cell image (ontresearch/wf-single-cell on Docker Hub), which runs the single-cell
# steps of popgen-lr-single-cell, ONT's wf-single-cell under Nextflow on Seqera. Multi-tool and built
# by ONT, so it is copied, not rebuilt. The tag is the SHA upstream pins as params.wf.container_sha,
# kept as is so the pipeline can build the image path from it.
SOURCE_IMAGE="docker://docker.io/ontresearch/wf-single-cell:sha69b55755e85ed2abf209bf19d9f72d7b5dd5f807"
IMAGE_NAME="wf-single-cell"
IMAGE_TAG="sha69b55755e85ed2abf209bf19d9f72d7b5dd5f807"

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy ${SOURCE_IMAGE} docker://australia-southeast1-docker.pkg.dev/cpg-common/images/${IMAGE_NAME}:${IMAGE_TAG}
