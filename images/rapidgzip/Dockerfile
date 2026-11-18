FROM python:3.12-slim-bookworm

# ENV, not ARG: the manual deploy workflows pass the image tag as the VERSION build argument, and a tag
# such as 0.16.2-1 is not a version pip can install.
ENV VERSION=0.16.2

LABEL description="rapidgzip: gzip decompression on many threads"

# rapidgzip is distributed on PyPI as a wheel that carries the command-line tool. procps lets Nextflow
# record task metrics. No USER is set, so the image runs as root: Hail Batch cannot run an image whose
# USER is a bare numeric UID.
RUN apt-get update && apt-get install -y --no-install-recommends procps && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    pip install --no-cache-dir rapidgzip==${VERSION}
