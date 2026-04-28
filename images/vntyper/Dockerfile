FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.137.cpg1-2 AS base

ARG VERSION=2.0.3

ENV MAMBA_ROOT_PREFIX=/root/micromamba
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH

# Install micromamba and the main vntyper dependencies
RUN apt-get update && apt-get install -y wget bzip2 && \
    rm -rf /var/lib/apt/lists/* && \
    wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/2.5.0/linux-64/micromamba-2.5.0-1.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir -p ${MAMBA_ROOT_PREFIX} && \
    micromamba install -y --prefix ${MAMBA_ROOT_PREFIX} -c bioconda -c conda-forge \
        bwa=0.7.18 \
        samtools=1.20 \
        bcftools=1.21 \
        fastp=0.23.4 && \
    rm -rf ${MAMBA_ROOT_PREFIX}/pkgs

# Create the envadvntr conda environment for adVNTR (Python 2.7)
COPY vntyper/conda/environment_envadvntr.yml /tmp/environment_envadvntr.yml
RUN micromamba create -y -f /tmp/environment_envadvntr.yml && \
    rm -rf ${MAMBA_ROOT_PREFIX}/pkgs /tmp/environment_envadvntr.yml

# Install adVNTR into the envadvntr environment
RUN git clone https://github.com/berntpopp/adVNTR.git --branch enhanced_hmm /tmp/adVNTR && \
    cd /tmp/adVNTR && \
    micromamba run -n envadvntr \
        bash -c 'CFLAGS="${CFLAGS:-} -Wno-error=incompatible-pointer-types -Wno-error=int-conversion" python setup.py install' && \
    rm -rf /tmp/adVNTR

# Install vntyper from source and download hg38 + advntr references
RUN git clone https://github.com/hassansaei/vntyper.git /tmp/vntyper-src && \
    cd /tmp/vntyper-src && \
    git checkout v${VERSION} && \
    pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir . && \
    vntyper install-references --output-dir /vntyper/reference --references hg38 advntr && \
    rm -rf /tmp/vntyper-src

# Patch vntyper to use the CRAM_REFERENCE environment variable
RUN sed -i 's|cram_ref_option = ""|cram_ref_option = f"--reference {os.environ['"'"'CRAM_REFERENCE'"'"']}" if os.environ.get("CRAM_REFERENCE") else ""|' \
    /usr/local/lib/python3.11/site-packages/vntyper/scripts/fastq_bam_processing.py

# Create symlinks so vntyper's relative paths resolve
RUN ln -s /usr/local/lib/python3.11/site-packages/vntyper /vntyper/vntyper && \
    ln -s ${MAMBA_ROOT_PREFIX}/bin/micromamba /usr/local/bin/mamba

WORKDIR /vntyper
