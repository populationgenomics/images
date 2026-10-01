FROM debian:bookworm-slim AS builder

# 1.20 is the last IBIS release; this pins the ahtremblay-patch-recommendation
# branch tip, 30 commits on, which fixes a variable type cast in the class-file
# output. The genetio submodule is pinned by that commit.
ENV VERSION=1.20-0e1d9aa
ENV IBIS_COMMIT=0e1d9aa3057d40b018feb0da488420901ba40f58

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        g++ \
        git \
        make \
        zlib1g-dev && \
    git clone https://github.com/williamslab/ibis.git /ibis && \
    cd /ibis && \
    git checkout "${IBIS_COMMIT}" && \
    git submodule update --init && \
    make -j"$(nproc)"

FROM debian:bookworm-slim

ENV VERSION=1.20-0e1d9aa

LABEL description="IBIS IBD segment caller for unphased array data, with seg2coef and bseg2seg"

# libgomp1: IBIS is built with -fopenmp
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgomp1 \
        zlib1g && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

COPY --from=builder /ibis/ibis /ibis/seg2coef /ibis/bseg2seg /ibis/maplen.awk /usr/local/bin/
COPY --from=builder /ibis/LICENSE /usr/local/share/doc/ibis/LICENSE
