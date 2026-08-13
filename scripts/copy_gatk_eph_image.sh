#!/usr/bin/env bash
set -ex
# Experimental GATK build for the FederateSVs tool, suggested by
# Emma Pierce-Hoffman (Talkowski lab).
# Upstream: https://github.com/broadinstitute/gatk
#   branch eph_federate_svs @ 53b4f16a7fd07cc66b7560e6beaaa7fbe80d61ca
# NOTE: the source lives in Broad's ephemeral (eph) registry and will be
# deleted eventually. Rebuild from the commit above with a Dockerfile
# before using this for any production callset.
SOURCE_IMAGE="docker://us.gcr.io/broad-dsde-methods/eph/gatk:eph_federate_svs-53b4f16a7"
IMAGE_NAME="gatk_federate_svs"
IMAGE_TAG="53b4f16a7"

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy ${SOURCE_IMAGE} docker://australia-southeast1-docker.pkg.dev/tenk10k-sv/images/${IMAGE_NAME}:${IMAGE_TAG}
