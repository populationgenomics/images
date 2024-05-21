FROM debian:buster-slim

# Set environment variables for Micromamba
ENV MAMBA_ROOT_PREFIX /root/micromamba
ENV PATH $MAMBA_ROOT_PREFIX/bin:$PATH

# Specify version and Git repository details as build arguments
ARG VERSION=${VERSION:-6.0.1}
ARG GIT_REPO=https://github.com/populationgenomics/TRTools.git
ARG GIT_BRANCH=skip-if-monomorphic

# Install necessary packages and clean up to reduce image size
RUN apt-get update && apt-get install -y \
    git wget bash bzip2 zip && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/cache/apt/* && \
    wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir ${MAMBA_ROOT_PREFIX}

# Clone the specified branch from the GitHub repository
RUN git clone -b ${GIT_BRANCH} ${GIT_REPO}

# Set the working directory
WORKDIR TRTools

# Create the Conda environment and install dependencies using Poetry
# COPY dev-env.yml .
RUN micromamba create -n trtools -f dev-env.yml && \
    micromamba run -n trtools poetry install

# Activate the Conda environment
CMD ["micromamba", "run", "-n", "trtools", "poetry", "shell"]
