FROM ubuntu:22.04

ARG VERSION=${VERSION:-14.0.0}

# download & unzip the distribution - this will create a directory called 'exomiser-cli-14.0.0'
RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    openjdk-21-jre-headless \
    python3 \
    python3-pip \
    tree \
    unzip \
    wget \
    && \
    apt-get clean; \
    wget "https://github.com/exomiser/Exomiser/releases/download/${VERSION}/exomiser-cli-${VERSION}-distribution.zip" &&  \
    unzip "exomiser-cli-${VERSION}-distribution.zip" && \
    rm "exomiser-cli-${VERSION}-distribution.zip" && \
    pip3 install PyYAML

WORKDIR /exomiser-cli-${VERSION}

COPY application.properties /exomiser-cli-${VERSION}/application.properties
