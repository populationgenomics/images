#!/usr/bin/env bash
set -ex
# ONT's shared tools image (ontresearch/wf-common on Docker Hub), which runs ingress (fastcat), BAM
# merging and the report in popgen-lr-single-cell, ONT's wf-single-cell under Nextflow on Seqera.
# Multi-tool and built by ONT, so it is copied, not rebuilt. The tag is the SHA upstream pins as
# params.wf.common_sha, kept as is so the pipeline can build the image path from it.
SOURCE_IMAGE="docker://docker.io/ontresearch/wf-common:shafdd79f8e4a6faad77513c36f623693977b92b08e"
IMAGE_NAME="wf-common"
IMAGE_TAG="shafdd79f8e4a6faad77513c36f623693977b92b08e"

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy ${SOURCE_IMAGE} docker://australia-southeast1-docker.pkg.dev/cpg-common/images/${IMAGE_NAME}:${IMAGE_TAG}
