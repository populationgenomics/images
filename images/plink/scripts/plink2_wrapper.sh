#!/bin/sh

# PLINK2 Multi-Architecture Wrapper
# Selects the best binary variant based on host CPU features

# Default to generic
BINARY="/usr/local/bin/variants/plink2_generic"

if [ -f /proc/cpuinfo ]; then
    # Check for AVX2 support
    if grep -q "avx2" /proc/cpuinfo; then
        # Check for AMD vs Intel optimization
        if grep -q "AuthenticAMD" /proc/cpuinfo; then
            BINARY="/usr/local/bin/variants/plink2_amd_avx2"
        else
            BINARY="/usr/local/bin/variants/plink2_intel_avx2"
        fi
    fi
fi

# Log the selection for debugging
# echo "Launching PLINK2 variant: $BINARY" >&2

exec "$BINARY" "$@"
