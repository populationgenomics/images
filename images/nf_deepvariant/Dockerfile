ARG VERSION=1.10.0

FROM google/deepvariant:${VERSION}

ARG VERSION

LABEL \
  version="v${VERSION}" \
  description="DeepVariant with procps for NextFlow"

# procps provides ps, which NextFlow needs to record task CPU and memory use
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*
