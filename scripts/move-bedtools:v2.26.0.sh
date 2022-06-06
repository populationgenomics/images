#!/usr/bin/env bash

set -ex

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy docker://quay.io/biocontainers/bedtools docker://australia-southeast1-docker.pkg.dev/cpg-common/images/bedtools:v2.26.0
