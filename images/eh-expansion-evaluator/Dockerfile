FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

ARG VERSION=${VERSION:-eb982ae41f83aaf44d4a8142980f7a41ac33e979}

# Install Git to clone the repository
RUN apt-get update && \
    apt-get install -y \
    git \
    && \
    apt-get clean -y

# Clone the repository
RUN git clone https://github.com/oscarlr/EH-expansion-evaluator

# Container should start in git repository
WORKDIR EH-expansion-evaluator
