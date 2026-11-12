# SHAPEIT5 (phase_common, phase_rare, ligate, switch, simulate, xcftools), built from source at a
# pinned commit of https://github.com/odelaneau/shapeit.
#
# Why this commit: the original odelaneau/shapeit5 repository is blocked by GitHub and its
# v5.1.1 release binaries are gone. The author's replacement repository, odelaneau/shapeit, has
# no releases or tags, only the single commit c34d4db3e99a2f7e23deb727671ae260901a5886
# (2026-05-19). That commit is newer than the v5.1.1 release (it adds xcftools and other changes),
# although versions/versions.h still reports "1.1" for each phasing tool, so the tools' version
# banners print "5.1.1" followed by "commit = c34d4db". The image tag is therefore the short
# commit, not a SHAPEIT5 release number.
FROM debian:trixie-slim AS base

# CI's get_version.py derives the image tag by text-matching the first
# VERSION assignment in this file (comments included), so keep the value
# below literally spelt out, with no variables. It must be the first 7 characters
# of SHAPEIT5_COMMIT; checked in the builder stage.
ENV VERSION=c34d4db
ENV SHAPEIT5_COMMIT=c34d4db3e99a2f7e23deb727671ae260901a5886
LABEL maintainer="Population Genomics Team"
LABEL version="${VERSION}"
LABEL org.opencontainers.image.source="https://github.com/odelaneau/shapeit"
LABEL org.opencontainers.image.revision="${SHAPEIT5_COMMIT}"

# Build stage: the compiler, libraries and sources stay out of the final image.
FROM base AS builder

# Must match the date of SHAPEIT5_COMMIT above; compiled into each tool's version banner.
ARG SHAPEIT5_COMMIT_DATE=2026-05-19
# sha256 of the GitHub archive tarball of SHAPEIT5_COMMIT (computed once with
# `curl -L <url> | sha256sum`), so the build fails if the archive changes.
ARG SHAPEIT5_SHA256=dcdad851e86ba01d8ddbf330e426f97c16e40cf64950dfd25d0641f1cd5f61f8
# htslib is built from source because Debian's libhts.a is not usable for a static link:
# libhts-dev depends on libcurl4-gnutls-dev, whose static libraries need the whole GnuTLS
# stack. The author's xcftools makefile references htslib-1.21 and the upstream Dockerfile
# uses 1.20.
ARG HTSLIB_VERSION=1.21
ARG HTSLIB_SHA256=84b510e735f4963641f26fd88c8abdee81ff4cb62168310ae716636aac0f1823

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bzip2 \
        ca-certificates \
        g++ \
        gcc \
        libboost-dev \
        libboost-iostreams-dev \
        libboost-program-options-dev \
        libbz2-dev \
        libdeflate-dev \
        liblzma-dev \
        make \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# Static htslib without libcurl, S3, GCS or plugins, so the SHAPEIT5 tools can read local
# files only (not http(s), ftp, gs:// or s3:// URLs); Hail Batch localises inputs before the
# tools run. No CPU flags are set here: htslib's SIMD code is chosen at run time.
RUN set -eu; \
    wget --quiet "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"; \
    printf '%s  %s\n' "${HTSLIB_SHA256}" "htslib-${HTSLIB_VERSION}.tar.bz2" > htslib.sha256; \
    sha256sum --check --strict htslib.sha256; \
    tar -xjf "htslib-${HTSLIB_VERSION}.tar.bz2"; \
    cd "htslib-${HTSLIB_VERSION}"; \
    ./configure --disable-libcurl --disable-s3 --disable-gcs --disable-plugins; \
    make -j"$(nproc)"; \
    make install; \
    mkdir -p /install/doc; \
    install -m 0644 LICENSE /install/doc/htslib-LICENSE

# Skip the repository's test data (about 360 MB), reference resources and website; none is
# needed to build.
RUN set -eu; \
    case "${SHAPEIT5_COMMIT}" in \
        "${VERSION}"*) ;; \
        *) echo "SHAPEIT5_COMMIT=${SHAPEIT5_COMMIT} does not start with VERSION=${VERSION}" >&2; exit 1 ;; \
    esac; \
    wget --quiet --output-document=shapeit.tar.gz \
        "https://github.com/odelaneau/shapeit/archive/${SHAPEIT5_COMMIT}.tar.gz"; \
    printf '%s  %s\n' "${SHAPEIT5_SHA256}" shapeit.tar.gz > shapeit.sha256; \
    sha256sum --check --strict shapeit.sha256; \
    mkdir shapeit; \
    tar -xzf shapeit.tar.gz --strip-components=1 -C shapeit \
        --exclude="shapeit-${SHAPEIT5_COMMIT}/test" \
        --exclude="shapeit-${SHAPEIT5_COMMIT}/resources" \
        --exclude="shapeit-${SHAPEIT5_COMMIT}/docs"; \
    rm shapeit.tar.gz; \
    install -m 0644 shapeit/LICENSE /install/doc/shapeit5-LICENSE

# Upstream's own `static_exe` make target (the one behind its Docker image and static
# release binaries) links the tools statically: libhts, Boost, zlib, bzip2, lzma and
# libdeflate go into the executable, so the final image needs no runtime libraries and the
# binaries do not depend on the image's libc.
#
# Overridden on the make command line rather than patching the sources:
#  - HTSLIB_INC / HTSLIB_LIB: the makefile points at a htslib-1.18 directory under an unset
#    HTSSRC; use the htslib built above in /usr/local.
#  - BOOST_LIB_IO / BOOST_LIB_PO: Debian's static Boost libraries. The common makefile of
#    the five phasing tools already points at these paths for static_exe, but xcftools'
#    makefile instead searches with `whereis`, which does not look in Debian's multiarch
#    directory, and the link then fails with undefined Boost references.
#  - DYN_LIBS_FOR_STATIC: drop -lcurl and -lcrypto (htslib is built without libcurl).
#  - COMMIT_VERS / COMMIT_DATE: the makefile runs `git`, which finds no repository in an
#    archive tarball and would leave the commit blank in the version banner.
#
# CPU flags: the static_exe target hardcodes -O2 -mavx2 -mfma (the other targets use -mavx2
# -mfma too; none uses -march=native). The binaries therefore need an AVX2- and FMA-capable
# x86-64 CPU (Haswell or later, AMD Excavator or later) and do not depend on the CPU of the
# machine that builds the image. Older n1 platforms (Sandy Bridge, Ivy Bridge) lack AVX2, so
# a job on a VM without it dies with SIGILL rather than an error message.
#
# Each tool is installed under its plain name, plus the <tool>_static name that upstream's
# documentation and Docker image use.
RUN set -eu; \
    mkdir -p /install/bin; \
    for tool in phase_common phase_rare ligate switch simulate xcftools; do \
        make -C "shapeit/${tool}" -j"$(nproc)" static_exe \
            COMMIT_VERS="${VERSION}" \
            COMMIT_DATE="${SHAPEIT5_COMMIT_DATE}" \
            HTSLIB_INC=/usr/local/include \
            HTSLIB_LIB=/usr/local/lib/libhts.a \
            BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a \
            BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a \
            DYN_LIBS_FOR_STATIC="-lz -lpthread -lbz2 -llzma -ldeflate"; \
        install -m 0755 -s "shapeit/${tool}/bin/${tool}_static" "/install/bin/${tool}"; \
        ln -s "${tool}" "/install/bin/${tool}_static"; \
    done

FROM base
COPY --from=builder /install/bin/ /usr/local/bin/
COPY --from=builder /install/doc/ /usr/share/doc/shapeit5/
