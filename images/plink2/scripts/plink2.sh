#!/bin/sh
# Run plink2, ldsc, mtag or vcf_subset, using the build that matches the host CPU.
#
# Installed as /usr/local/bin/plink2, with ldsc, mtag and vcf_subset symlinked to it;
# the tool is taken from the name this script was called as. Each tool is installed as
# /usr/local/bin/variants/<tool>_<variant> for the variants generic, intel_avx2 and
# amd_avx2.
#
# PLINK2_FORCE_VARIANT (generic | intel_avx2 | amd_avx2) overrides CPU detection,
# e.g. to sidestep an AVX2-build-specific bug.

set -eu

TOOL=$(basename "$0")
VARIANT_DIR=/usr/local/bin/variants

if [ -n "${PLINK2_FORCE_VARIANT:-}" ]; then
    BINARY="$VARIANT_DIR/${TOOL}_$PLINK2_FORCE_VARIANT"
    if [ ! -x "$BINARY" ]; then
        echo "PLINK2_FORCE_VARIANT=$PLINK2_FORCE_VARIANT: $BINARY not found" >&2
        exit 1
    fi
else
    VARIANT=generic
    if grep -q "avx2" /proc/cpuinfo 2>/dev/null; then
        if grep -q "AuthenticAMD" /proc/cpuinfo; then
            VARIANT=amd_avx2
        else
            VARIANT=intel_avx2
        fi
    fi
    BINARY="$VARIANT_DIR/${TOOL}_$VARIANT"
    if [ ! -x "$BINARY" ]; then
        echo "$TOOL: $BINARY not found" >&2
        exit 1
    fi
fi

exec "$BINARY" "$@"
