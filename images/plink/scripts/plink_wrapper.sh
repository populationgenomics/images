#!/bin/sh

# PLINK Multi-Version Wrapper
# Handles plink (1.9), plink2 (2.0), and python commands

# Get the base name of the command used to call this script
BASE=$(basename "$0")
CMD="$1"

# Select the best PLINK2 binary variant based on host CPU features
select_plink2_binary() {
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
    echo "$BINARY"
}

# Routing logic
# 1. If called as 'plink2' OR first arg is 'plink2'
if [ "$BASE" = "plink2" ] || [ "$CMD" = "plink2" ]; then
    if [ "$CMD" = "plink2" ]; then shift; fi
    BINARY=$(select_plink2_binary)
    # echo "[WRAPPER] Routing to PLINK2: $BINARY" >&2
    exec "$BINARY" "$@"

# 2. If called as 'plink' OR first arg is 'plink' or 'plink1.9'
elif [ "$BASE" = "plink" ] || [ "$CMD" = "plink" ] || [ "$CMD" = "plink1.9" ]; then
    if [ "$CMD" = "plink" ] || [ "$CMD" = "plink1.9" ]; then shift; fi
    # echo "[WRAPPER] Routing to PLINK 1.9" >&2
    exec /usr/local/bin/plink1.9 "$@"

# 3. If first arg is python/python3
elif [ "$CMD" = "python" ] || [ "$CMD" = "python3" ]; then
    # echo "[WRAPPER] Routing to Python" >&2
    exec "$@"

# 4. Default / Fallback
else
    # If it starts with --, assume plink 1.9
    if echo "$CMD" | grep -q "^-"; then
        # echo "[WRAPPER] Defaulting to PLINK 1.9 (arg starts with --)" >&2
        exec /usr/local/bin/plink1.9 "$@"
    else
        # Try to execute as-is
        # echo "[WRAPPER] Executing as-is: $@" >&2
        exec "$@"
    fi
fi
