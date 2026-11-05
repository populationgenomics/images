FROM debian:trixie-slim

ARG VERSION=0.3.15

RUN apt-get update \
    && apt-get install -y --no-install-recommends  \
      procps \
    && rm -rf /var/lib/apt/lists/*

ADD https://github.com/brentp/mosdepth/releases/download/v${VERSION}/mosdepth  /usr/local/bin/

RUN chmod a+x /usr/local/bin/mosdepth
