ARG SVAFOTATE_COMMIT=30b5004a0f4d26959c6b9a82f165651585293626

ARG VERSION=0.0.1

FROM python:3.11-slim-bullseye AS base

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        bzip2 \
        ca-certificates \
        git \
        gnupg \
        libbz2-1.0 \
        libcurl4 \
        liblzma5 \
        procps \
        wget \
        zip \
        zlib1g && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

FROM base AS compiler

ARG BCFTOOLS_VERSION=1.24

RUN apt-get update && apt-get install --no-install-recommends -y \
        gcc \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libssl-dev \
        make \
        zlib1g-dev && \
    rm -rf /var/lib/apt/lists/* && \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -xf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    ./configure --enable-libcurl --enable-s3 --enable-gcs && \
    make && \
    strip bcftools plugins/*.so && \
    make DESTDIR=/bcftools_install install && \
    cd htslib-${BCFTOOLS_VERSION} && \
    make && \
    make DESTDIR=/bcftools_install install

# SVAFotate's dependency set is not wheel-complete on every architecture - sorted_nearest (via
# pyranges) publishes no wheels at all, and ncls publishes none for manylinux aarch64, so an arm64
# build compiles both from sdist. Wheeling the whole set here rather than naming the sdist-only
# packages keeps that architecture-dependent, and keeps gcc out of the final image.
RUN pip wheel --no-cache-dir --wheel-dir /sv_wheels \
        'pandas<3' \
        pyranges \
        ncls

FROM base AS talos

ARG SVAFOTATE_COMMIT

COPY --from=compiler /bcftools_install/usr/local/bin/* /usr/local/bin/
COPY --from=compiler /bcftools_install/usr/local/libexec/bcftools/* /usr/local/libexec/bcftools/

ARG ECHTVAR_VERSION=v0.2.2
ENV VERSION=0.0.1

RUN wget -q -O /bin/echtvar "https://github.com/brentp/echtvar/releases/download/${ECHTVAR_VERSION}/echtvar" && \
    chmod +x /bin/echtvar

COPY --from=ghcr.io/astral-sh/uv:0.9.26 /uv /uvx /bin/

# Enable bytecode compilation
ENV UV_COMPILE_BYTECODE=1

# Copy from the cache instead of linking since it's a mounted volume
ENV UV_LINK_MODE=copy

RUN git clone https://github.com/populationgenomics/talos2.git

WORKDIR /talos2

# Install the project's dependencies using the lockfile and settings
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --all-extras --frozen --no-install-project --no-dev

COPY --from=compiler /sv_wheels /tmp/sv_wheels
RUN --mount=type=cache,target=/root/.cache/uv \
    uv pip install --python /talos/.venv/bin/python --find-links /tmp/sv_wheels \
        'pandas<3' \
        pyranges \
        ncls \
        "https://github.com/fakedrtom/SVAFotate/archive/${SVAFOTATE_COMMIT}.tar.gz" && \
    rm -rf /tmp/sv_wheels

# Place executables in the environment at the front of the path
ENV PATH="/talos2/.venv/bin:$PATH"

RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --all-extras --frozen --no-dev --inexact

# both halves of the image have to survive each other's install - the Talos package and its cyvcf2,
# and the SVAFotate CLI. Cheap, and the failure mode it catches is otherwise a runtime one
RUN python -c "import talos2, cyvcf2, svafotate" && \
    bcftools --version && \
    tabix --version && \
    echtvar --version
