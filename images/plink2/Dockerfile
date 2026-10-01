# PLINK 2.0 development build (beta 1, 20260929) with ldsc, mtag, vcf_subset, bgzip and tabix.
#
# This is the dev build, not the stable alpha release in the plink image: it lifts the
# restriction that --pmerge and --pmerge-list could only do concatenation-like jobs.
# Its downloads therefore come from the plink2-assets root, NOT the alpha7/ subdirectory.
FROM debian:trixie-slim AS base

# CI's get_version.py derives the image tag by text-matching the first
# VERSION assignment in this file (comments included), so keep the value
# below literally spelt out, with no variables.
ENV VERSION=2.0.0-b.1-20260929
LABEL maintainer="Population Genomics Team"
LABEL version="${VERSION}"
# ldsc and mtag are plink-ng's own v0.1 C reimplementations, NOT the reference Python
# LDSC (Bulik-Sullivan) or MTAG (Turley) tools; their options and output can differ.
LABEL ldsc_mtag_note="ldsc and mtag are plink-ng's own v0.1 C reimplementations, not the reference Python LDSC (Bulik-Sullivan) or MTAG (Turley) tools"

# The PLINK binaries are statically linked; tabix provides tabix and bgzip.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        tabix \
    && rm -rf /var/lib/apt/lists/*

# Download stage: wget and unzip stay out of the final image.
FROM base AS downloader

# Must match the date at the end of VERSION above; checked below.
ARG PLINK2_BUILD=20260929
# sha256 of each zip, so the build fails if upstream replaces a file.
ARG SHA256_INTEL_AVX2=2b667c308538c7e59e37c393d79b945bf58181c7e7bcf1f0f1bf989e9071ee9b
ARG SHA256_AMD_AVX2=c07b77f87998ed887ebabf76d3f41070ad078af285607ca9508f1f1ece362d03
ARG SHA256_GENERIC=6149fe9477d7a5cbbbc74f1e9c5c972939a6ade313484ff12c13c88e610652e5

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        unzip \
        wget \
    && rm -rf /var/lib/apt/lists/*

# Each zip holds plink2, ldsc, mtag and vcf_subset. The tools are built per CPU variant
# (the intel/amd AVX2 ldsc, mtag and vcf_subset are byte-identical to each other but
# differ from the generic ones), so all four are installed per variant as
# /install/variants/<tool>_<variant> and chosen at run time by scripts/plink2.sh.
# Only the Intel and generic zips contain the Intel license file; they are identical.
WORKDIR /tmp/plink2
RUN set -eu; \
    case "${VERSION}" in \
        *-"${PLINK2_BUILD}") ;; \
        *) echo "VERSION=${VERSION} does not end in PLINK2_BUILD=${PLINK2_BUILD}" >&2; exit 1 ;; \
    esac; \
    mkdir -p /install/variants /install/doc; \
    install_variant() { \
        zip_variant="$1"; variant="$2"; sha256="$3"; \
        zip="plink2_linux_${zip_variant}_${PLINK2_BUILD}.zip"; \
        wget --quiet "https://s3.amazonaws.com/plink2-assets/${zip}"; \
        printf '%s  %s\n' "${sha256}" "${zip}" > "${zip}.sha256"; \
        sha256sum --check --strict "${zip}.sha256"; \
        unzip -q "${zip}" -d "${variant}"; \
        for tool in plink2 ldsc mtag vcf_subset; do \
            install -m 0755 "${variant}/${tool}" "/install/variants/${tool}_${variant}"; \
        done; \
    }; \
    install_variant avx2 intel_avx2 "${SHA256_INTEL_AVX2}"; \
    install_variant amd_avx2 amd_avx2 "${SHA256_AMD_AVX2}"; \
    install_variant x86_64 generic "${SHA256_GENERIC}"; \
    install -D -m 0644 intel_avx2/intel-simplified-software-license.txt \
        /install/doc/intel-simplified-software-license.txt

FROM base
COPY --from=downloader /install/variants/ /usr/local/bin/variants/
COPY --from=downloader /install/doc/ /usr/share/doc/plink2/

# One dispatch script, installed as plink2 and symlinked as ldsc, mtag and vcf_subset;
# it picks the CPU variant of whichever tool it was called as.
COPY scripts/plink2.sh /usr/local/bin/plink2
RUN chmod 0755 /usr/local/bin/plink2 \
    && ln -s plink2 /usr/local/bin/ldsc \
    && ln -s plink2 /usr/local/bin/mtag \
    && ln -s plink2 /usr/local/bin/vcf_subset
