FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.138.cpg2-2

# A temporary image for testing the metamist 8.0.0 pre-releases in a Hail Batch pipeline:
# cpg_hail_gcloud plus a pinned metamist client.
# ENV rather than ARG: the base image sets ENV VERSION (which would override an ARG of
# the same name), and this keeps the pin immune to the CI `--build-arg VERSION=<tag>-<n>`.
ENV VERSION=8.0.0.dev65

RUN pip install --no-cache-dir "metamist==${VERSION}" && \
    python -c "import hail, metamist" && \
    gcloud --version
