FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.137.cpg1-2 AS base

ARG VERSION=${VERSION:-2.0.3}

ENV MAMBA_ROOT_PREFIX=/root/micromamba
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH

# Install micromamba and use it to install vntyper dependencies
RUN apt-get update && apt-get install -y wget bzip2 && \
    rm -rf /var/lib/apt/lists/* && \
    wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir -p ${MAMBA_ROOT_PREFIX} && \
    micromamba install -y --prefix ${MAMBA_ROOT_PREFIX} -c bioconda -c conda-forge \
        bwa=0.7.18 \
        samtools=1.20 \
        bcftools=1.21 \
        fastp=0.23.4 && \
    rm -rf ${MAMBA_ROOT_PREFIX}/pkgs

# Install vntyper from source and download hg38 references
RUN git clone https://github.com/hassansaei/vntyper.git && \
    cd vntyper && \
    git checkout v${VERSION} && \
    pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir . && \
    vntyper install-references --output-dir /vntyper/reference --references hg38

# Patch vntyper to use the CRAM_REFERENCE environment variable, necessary for compatibility with CPG's CRAM reference management
RUN sed -i 's|cram_ref_option = ""|cram_ref_option = f"--reference {os.environ['"'"'CRAM_REFERENCE'"'"']}" if os.environ.get("CRAM_REFERENCE") else ""|' \
      /usr/local/lib/python3.11/site-packages/vntyper/scripts/fastq_bam_processing.py

# Create symlink for vntyper package install directory
RUN ln -s /usr/local/lib/python3.11/site-packages/vntyper /vntyper/vntyper

WORKDIR /vntyper
