FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.138.cpg2-1 AS base

ARG VERSION=${VERSION:-1.34}

RUN pip install multiqc==${VERSION}
