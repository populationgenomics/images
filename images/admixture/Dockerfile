# ADMIXTURE 1.4.0, the upstream static Linux binary, for maximum-likelihood ancestry estimation.
#
# Usage: admixture [-jN] [--cv] [--seed=X] <prefix>.bed <K>
# The <prefix>.K.Q and <prefix>.K.P outputs are written to the working directory, not next to
# the input.
FROM debian:trixie-slim AS downloader

# CI's get_version.py derives the image tag by text-matching the first
# VERSION assignment in this file (comments included), so keep the value
# below literally spelt out, with no variables. A version bump must also
# update the VERSION in the final stage and SHA256 below.
ENV VERSION=1.4.0
# sha256 of the release tarball, so the build fails if upstream replaces it.
ARG SHA256=5520c4f18525f49ab2c715d711b40f4a0818ef1485361a94ec8954d1364e6f19

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        wget && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /tmp/admixture
RUN wget --quiet "https://dalexander.github.io/admixture/binaries/admixture_linux-${VERSION}.tar.gz" && \
    printf '%s  %s\n' "${SHA256}" "admixture_linux-${VERSION}.tar.gz" | sha256sum --check --strict && \
    tar -xzf "admixture_linux-${VERSION}.tar.gz"

FROM debian:trixie-slim

ENV VERSION=1.4.0
LABEL maintainer="Population Genomics Team"
LABEL version="${VERSION}"
LABEL description="ADMIXTURE maximum-likelihood ancestry estimation from PLINK genotype files"

# The admixture binary is statically linked, so it needs no runtime libraries.
COPY --from=downloader /tmp/admixture/admixture_linux-${VERSION}/admixture /usr/local/bin/admixture
COPY --from=downloader /tmp/admixture/admixture_linux-${VERSION}/admixture-manual.pdf /usr/local/share/doc/admixture/admixture-manual.pdf

# Fail the build if the binary cannot execute; --help exits 0 (no arguments exits 255).
RUN admixture --help > /dev/null
