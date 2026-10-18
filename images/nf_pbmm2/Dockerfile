FROM australia-southeast1-docker.pkg.dev/cpg-common/images/nf_base:5.8.1-1

ARG VERSION=26.2.0

LABEL \
  version="v${VERSION}" \
  description="NF Base pbmm2 install"

ENV MAMBA_ROOT_PREFIX=/opt/micromamba

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# pbmm2 lives in its own conda env; only its binary goes on PATH, so nf_base's python stays first
RUN wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/2.9.0/linux-64/micromamba-2.9.0-0.tar.bz2 | \
        tar -xj -C /usr/local bin/micromamba && \
    micromamba create -y --prefix /opt/pbmm2 -c conda-forge -c bioconda pbmm2=${VERSION} && \
    micromamba clean -y --all && \
    ln -s /opt/pbmm2/bin/pbmm2 /usr/local/bin/pbmm2
