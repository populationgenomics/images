#!/usr/bin/env bash

set -ex

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy docker://quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0 docker://australia-southeast1-docker.pkg.dev/cpg-common/images/bedtools:v2.30.0
