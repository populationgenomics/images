FROM debian:bookworm-slim

ENV VERSION=0.1.1

LABEL description="Genotype-based single-cell demultiplexing: cellsnp-lite, vireoSNP, pysam, samtools/bcftools/bedtools, plus git, gcloud and cpg-utils for analysis-runner jobs"

ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV PATH=${MAMBA_ROOT_PREFIX}/bin:${PATH}

# Login shells re-read /etc/profile and would otherwise drop the ENV PATH.
RUN echo 'PATH=/opt/micromamba/bin:$PATH' > /etc/profile.d/micromamba.sh

# git: prepare_git_job clones the calling repo into the job.
# gcloud: jobs read and write dataset buckets directly.
RUN apt-get update && apt-get install -y --no-install-recommends \
        bzip2 \
        ca-certificates \
        curl \
        git \
        gnupg && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
        > /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
        gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    apt-get update && \
    apt-get install -y --no-install-recommends google-cloud-cli && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

# Pinned so a rebuild of this tag reproduces the same toolchain.
RUN curl -fsSL -o /usr/local/bin/micromamba \
        https://github.com/mamba-org/micromamba-releases/releases/download/1.5.8-0/micromamba-linux-64 && \
    chmod +x /usr/local/bin/micromamba && \
    micromamba install -y --prefix ${MAMBA_ROOT_PREFIX} \
        -c bioconda -c conda-forge \
        python=3.12 \
        cellsnp-lite=1.2.3 \
        vireosnp=0.5.9 \
        pysam=0.23.0 \
        samtools=1.21 \
        bcftools=1.21 \
        bedtools=2.31.1 \
        numpy \
        scipy && \
    micromamba clean -y --all && \
    rm -r ${MAMBA_ROOT_PREFIX}/pkgs

# prepare_git_job resolves the checkout-token secret by running config_retrieve
# inside the job. Its generated shell swallows import errors, so without
# cpg-utils here a private-repo clone fails with a misleading "no git
# credentials secret found".
RUN pip install --no-cache-dir cpg-utils
