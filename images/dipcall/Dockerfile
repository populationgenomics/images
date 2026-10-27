FROM debian:bookworm-slim

ARG VERSION=0.3

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# dipcall's release kit bundles its own minimap2, k8, samtools, htsbox and bedtk as linux x86-64
# binaries, plus hs38.PAR.bed. On top of it: perl for run-dipcall (Getopt::Std is in the full perl
# package, not perl-base), make to run the makefile it writes, tabix and bgzip (the makefile writes no
# VCF index and the caller renames the sample column), and procps, which Nextflow needs for task metrics.
RUN apt-get update && apt-get install -y --no-install-recommends \
        bzip2 \
        ca-certificates \
        make \
        perl \
        procps \
        tabix \
        wget && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    wget -qO- "https://github.com/lh3/dipcall/releases/download/v${VERSION}/dipcall-${VERSION}_x64-linux.tar.bz2" | \
        tar -xj -C /opt && \
    test -x /opt/dipcall.kit/run-dipcall

ENV PATH=/opt/dipcall.kit:$PATH

# Every binary the popgen-pangenome DIPCALL process calls resolves and runs.
RUN perl -c /opt/dipcall.kit/run-dipcall && \
    minimap2 --version && \
    samtools --version | head -n 1 && \
    k8 -v && \
    { htsbox 2>&1 || true; } | grep -q 'htsbox r' && \
    bedtk 2>&1 | head -n 1 && \
    tabix --version | head -n 1 && \
    bgzip --version | head -n 1 && \
    make --version | head -n 1 && \
    test -s /opt/dipcall.kit/hs38.PAR.bed
