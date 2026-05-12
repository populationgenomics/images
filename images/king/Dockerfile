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
RUN URL_VERSION=$(echo "${VERSION}" | tr -d '.') \
    && wget -q "https://www.kingrelatedness.com/executables/Linux-king${URL_VERSION}.tar.gz" -O /tmp/king.tar.gz \
    && tar -xzf /tmp/king.tar.gz -C /usr/local/bin king \
    && chmod +x /usr/local/bin/king \
    && rm /tmp/king.tar.gz
