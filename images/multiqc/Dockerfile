FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.134.cpg2-2 AS base

ARG VERSION=${VERSION:-1.32}

RUN pip install multiqc==${VERSION}
