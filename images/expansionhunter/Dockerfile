
FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG CONDA_RELEASE="4.6.14"
ARG CONDA_INSTALL_DIR="/opt/conda"
ARG CONDA_BIN=${CONDA_INSTALL_DIR}/bin
ARG CONDA_CMD=${CONDA_BIN}/conda
ENV PATH=${CONDA_BIN}:$PATH
ENV EH_VERSION=v5.0.0

# hadolint ignore=DL3008
RUN apt-get update && apt-get install --no-install-recommends -qqy \
    python3-dev=3.8.2-0ubuntu2 \
    python3-pip=20.0.2-5ubuntu1.6 \
    python3.9\
    python-dev \
    build-essential \
    git \
    libbz2-dev \
    liblzma-dev \
    make \
    pkg-config \
    wget \
    unzip \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev\
     && rm  -rf/ var/lib/apt/lists/*

RUN rm /var/lib/apt/lists/*

# Install and configure conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_RELEASE}-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p ${CONDA_INSTALL_DIR} && \
    ${CONDA_CMD} init bash && \
    bash -c "source ${HOME}/.bashrc"

RUN pip install --no-cache-dir Cython==0.29.24 numpy==1.21.2 pybedtools==0.8.2 scipy==1.7.1


# Install ExpansionHunter
RUN wget --progress=dot:giga https://github.com/Illumina/ExpansionHunter/releases/download/${EH_VERSION}/ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz && \
    tar xzf ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz && \
    rm ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz && \
    mv /ExpansionHunter-${EH_VERSION}-linux_x86_64 /ExpansionHunter
ENV PATH="/ExpansionHunter/bin/:$PATH"

# Clean up
RUN rm -rf /tmp/* \
    /var/tmp/* \
    /var/cache/apt/* \
    /var/lib/apt/lists/* \

