FROM debian:bookworm-slim

ARG VERSION=${VERSION:-0.3.11}

ADD https://github.com/brentp/mosdepth/releases/download/v${VERSION}/mosdepth  /usr/local/bin/

RUN chmod a+x /usr/local/bin/mosdepth