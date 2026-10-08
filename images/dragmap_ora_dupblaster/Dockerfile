FROM australia-southeast1-docker.pkg.dev/cpg-common/images/dragmap_plus_ora:2.7.0-2 AS build

# set bash as the shell of choice
SHELL ["/bin/bash", "-c"]

ARG VERSION=0.3.0

# install build-essential for compiling
# install Rust & Cargo to facilitate install of dupblaster
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        curl && \
    curl https://sh.rustup.rs -sSf | bash -s -- -y && \
    echo 'source $HOME/.cargo/env' >> $HOME/.bashrc && \
    source $HOME/.bashrc && \
    cargo install dupblaster --version ${VERSION}

FROM australia-southeast1-docker.pkg.dev/cpg-common/images/dragmap_plus_ora:2.7.0-2

COPY --from=build /root/.cargo/bin/dupblaster /usr/local/bin/dupblaster
