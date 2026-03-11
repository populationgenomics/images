# Multi-Version PLINK Image (1.9 & 2.0) with Python 3
FROM python:3.11-slim-bookworm

# Build metadata
ENV VERSION=1.9-20250819-PLINK-2.0-20260228
LABEL maintainer="Population Genomics Team"
LABEL version="${VERSION}"

# Install dependencies and cleanup in a single layer to minimize size
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    unzip \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Create directory for binary variants
RUN mkdir -p /usr/local/bin/variants

# 1. Download and install PLINK 1.9 (Build 20250819)
RUN wget -q https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip -O plink1.zip \
    && unzip -q plink1.zip -d /tmp/plink1 \
    && mv /tmp/plink1/plink /usr/local/bin/plink1.9 \
    && rm -rf plink1.zip /tmp/plink1

# 2. Download and install PLINK 2.0 variants (Build 20260228)
# Intel AVX2
RUN wget -q https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_20260228.zip -O intel.zip \
    && unzip -q intel.zip -d /tmp/intel \
    && mv /tmp/intel/plink2 /usr/local/bin/variants/plink2_intel_avx2 \
    && rm -rf intel.zip /tmp/intel

# AMD AVX2
RUN wget -q https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_amd_avx2_20260228.zip -O amd.zip \
    && unzip -q amd.zip -d /tmp/amd \
    && mv /tmp/amd/plink2 /usr/local/bin/variants/plink2_amd_avx2 \
    && rm -rf amd.zip /tmp/amd

# Generic x86_64
RUN wget -q https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20260228.zip -O generic.zip \
    && unzip -q generic.zip -d /tmp/generic \
    && mv /tmp/generic/plink2 /usr/local/bin/variants/plink2_generic \
    && rm -rf generic.zip /tmp/generic

# 3. Install the multi-version wrapper script
COPY src/popgen_genotyping/scripts/plink_wrapper.sh /usr/local/bin/plink-wrapper
RUN chmod +x /usr/local/bin/plink-wrapper /usr/local/bin/plink1.9 /usr/local/bin/variants/plink2_* && \
    ln -s /usr/local/bin/plink-wrapper /usr/local/bin/plink && \
    ln -s /usr/local/bin/plink-wrapper /usr/local/bin/plink2

# Set entrypoint to the wrapper
ENTRYPOINT ["/usr/local/bin/plink-wrapper"]
CMD ["plink1.9", "--help"]
