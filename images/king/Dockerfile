FROM python:3.11-slim-bookworm

ARG VERSION=2.3.2
ENV VERSION=${VERSION}
LABEL maintainer="Population Genomics Team"
LABEL version="${VERSION}"

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    libquadmath0 \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# KING upstream URL pattern strips dots from the version (e.g. 2.3.2 -> 232).
# Browser User-Agent + retries: kingrelatedness.com filters the default wget UA
# from some networks (notably CI runners), returning HTTP errors.
RUN URL_VERSION=$(echo "${VERSION}" | tr -d '.') \
    && wget --tries=3 --timeout=60 \
        --user-agent="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36" \
        "https://www.kingrelatedness.com/executables/Linux-king${URL_VERSION}.tar.gz" \
        -O /tmp/king.tar.gz \
    && tar -xzf /tmp/king.tar.gz -C /usr/local/bin king \
    && chmod +x /usr/local/bin/king \
    && rm /tmp/king.tar.gz
