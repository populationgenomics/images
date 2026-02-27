# ------------------------------------------------------------------------------
# Stage 1: Build Environment
# ------------------------------------------------------------------------------
# We use the official Rust image as our builder. It provides cargo, rustc, etc.
FROM rust:1.93-trixie AS builder

# htslib (and its underlying libz-sys) requires cmake, a C/C++ compiler,
# and zlib development headers to build successfully.
RUN apt-get update && apt-get install -y \
    cmake \
    build-essential \
    git \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libclang-dev \
    && rm -rf /var/lib/apt/lists/*

# git clone and compile the application in release mode
## This will output the compiled binary at target/release/dupmark
RUN git clone https://github.com/populationgenomics/rust_dupmark.git && \
    cd rust_dupmark && \
    cargo build --release

# ------------------------------------------------------------------------------
# Stage 2: Minimal Runtime Environment
# ------------------------------------------------------------------------------
# We use debian:trixie-slim to dramatically reduce the final image size.
# It strips out the Rust toolchain and unnecessary cache layers.
FROM debian:trixie-slim

# Install the runtime shared libraries that the compiled htslib may depend on
RUN apt-get update && apt-get install -y \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libcurl4 \
    && rm -rf /var/lib/apt/lists/*

# Copy the compiled executable from the builder stage
COPY --from=builder /rust_dupmark/target/release/dupmark /usr/local/bin/dupmark

RUN dupmark --help
