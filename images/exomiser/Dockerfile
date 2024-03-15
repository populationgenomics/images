FROM ubuntu:22.04

ARG VERSION=${VERSION:-14.0.0}

# download & unzip the distribution - this will create a directory called 'exomiser-cli-14.0.0'
RUN apt-get update && \
    apt-get install -y \
    openjdk-21-jre-headless \
    unzip \
    wget \
    && \
    apt-get clean; \
    wget "https://github.com/exomiser/Exomiser/releases/download/${VERSION}/exomiser-cli-${VERSION}-distribution.zip" &&  \
    unzip "exomiser-cli-${VERSION}-distribution.zip" && \
    rm "exomiser-cli-${VERSION}-distribution.zip"

WORKDIR /exomiser-cli-${VERSION}
