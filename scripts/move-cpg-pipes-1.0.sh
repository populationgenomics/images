#!/usr/bin/env bash

set -ex

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy docker://australia-southeast1-docker.pkg.dev/analysis-runner/images/cpg-pipes:1.0 docker://australia-southeast1-docker.pkg.dev/cpg-common/images/cpg-pipes:1.0
