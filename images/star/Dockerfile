FROM debian:bookworm-slim

ARG VERSION=${VERSION:-2.7.10b}

RUN apt-get update && apt-get install -y git wget bash bzip2 zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    cd /tmp && wget https://github.com/alexdobin/STAR/releases/download/${VERSION}/STAR_${VERSION}.zip && \
    unzip STAR_${VERSION}.zip && \
    mv STAR_${VERSION}/Linux_x86_64_static/STAR* /usr/local/bin/ && \
    chmod a+x /usr/local/bin/STAR* && \
    rm -rf /tmp/STAR_${VERSION}
