# Minimal Multi-Arch PLINK2 Image
FROM debian:stable-slim

# Build metadata
ENV VERSION=2.0-alpha-20260307
LABEL maintainer="Population Genomics Team"
LABEL version="${VERSION}"

# Install runtime and build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    unzip \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Create directory for binary variants
RUN mkdir -p /usr/local/bin/variants

# Download and install PLINK2 variants (Build ${VERSION})
# 1. Intel AVX2
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_20260307.zip -O intel.zip \
    && unzip intel.zip -d /tmp/intel \
    && mv /tmp/intel/plink2 /usr/local/bin/variants/plink2_intel_avx2 \
    && rm -rf intel.zip /tmp/intel

# 2. AMD AVX2
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_amd_avx2_20260307.zip -O amd.zip \
    && unzip amd.zip -d /tmp/amd \
    && mv /tmp/amd/plink2 /usr/local/bin/variants/plink2_amd_avx2 \
    && rm -rf amd.zip /tmp/amd

# 3. Generic x86_64
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20260307.zip -O generic.zip \
    && unzip generic.zip -d /tmp/generic \
    && mv /tmp/generic/plink2 /usr/local/bin/variants/plink2_generic \
    && rm -rf generic.zip /tmp/generic

# Install the wrapper script
COPY scripts/plink2_wrapper.sh /usr/local/bin/plink2
RUN chmod +x /usr/local/bin/plink2 /usr/local/bin/variants/plink2_*

# Set entrypoint to the wrapper
ENTRYPOINT ["/usr/local/bin/plink2"]
CMD ["--help"]
