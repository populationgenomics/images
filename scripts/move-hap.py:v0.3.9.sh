#!/usr/bin/env bash

set -ex

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy docker://docker.io/pkrusche/hap.py:v0.3.9 docker://australia-southeast1-docker.pkg.dev/cpg-common/images/hap.py:v0.3.9
