FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
ARG VERSION=${VERSION:-eb982ae41f83aaf44d4a8142980f7a41ac33e979}

# Install Git to clone the repository
RUN apt-get update && apt-get install --no-install-recommends -qqy \
    python3-dev=3.8.2-0ubuntu2 \
    python3-pip=20.0.2-5ubuntu1.9 \
    python3.6\
    git \
    && \
    apt-get clean -y  \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip3 install --no-cache-dir numpy==1.24.4 pandas==2.0.3 scikit-learn==1.3.0

# Clone the repository
RUN git clone https://github.com/oscarlr/EH-expansion-evaluator && cd EH-expansion-evaluator && git checkout $VERSION

# Container should start in git repository
WORKDIR EH-expansion-evaluator
