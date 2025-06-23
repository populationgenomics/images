# Use official Google DeepVariant image as base
ARG VERSION=1.9.0

# Choose CPU or GPU variant
FROM google/deepvariant:${VERSION} as deepvariant_base

# Final image with DeepVariant + additional tools
FROM deepvariant_base
LABEL maintainer="https://github.com/google/deepvariant/issues"