FROM australia-southeast1-docker.pkg.dev/cpg-common/images/dragmap:1.3.0-tokenizer-next-fix-1 as builder

# set bash as the shell of choice
SHELL ["/bin/bash", "-c"]

# install build-essential for compiling
# install Rust & Cargo to facilitate install of dupblaster
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        bzip2 \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        wget \
        zlib1g-dev && \
    curl https://sh.rustup.rs -sSf | bash -s -- -y && \
    echo 'source $HOME/.cargo/env' >> $HOME/.bashrc && \
    source $HOME/.bashrc && \
    cargo install dupblaster

FROM australia-southeast1-docker.pkg.dev/cpg-common/images/dragmap:1.3.0-tokenizer-next-fix-1 a

COPY --from=build /root/.cargo/bin/dupblaster /usr/local/bin/dupblaster

ARG VERSION=2.7.0

# source
# https://sapac.support.illumina.com/downloads/dragen_decompression_software-linux.html

# installation
# https://support-docs.illumina.com/SW/Decompression_v2_6/Content/SW/ORA/Decompression_Installation.htm

RUN wget -q https://s3.amazonaws.com/webdata.illumina.com/downloads/software/dragen-decompression/orad.2.7.0.linux.tar.gz && \
    tar -xzf orad.2.7.0.linux.tar.gz && \
    mv orad.2.7.0.linux/orad /usr/local/bin/orad && \
    rm -rf orad.2.7.0.linux.tar.gz orad.2.7.0.linux
