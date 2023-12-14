FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y no-install-recommends \
    apt-transport-https ca-certificates gnupg curl sudo

## Install Google cloud-sdk
RUN curl https://sdk.cloud.google.com > install.sh && \
    bash install.sh --disable-prompts --install-dir=/opt && \
    rm install.sh

ENV PATH=$PATH:/opt/google-cloud-sdk/bin
