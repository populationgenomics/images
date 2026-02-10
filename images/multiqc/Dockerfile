FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.137.cpg1-2 AS base

ARG VERSION=${VERSION:-1.33}

RUN pip install multiqc==${VERSION}
