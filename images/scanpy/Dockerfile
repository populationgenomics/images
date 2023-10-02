ARG DRIVER_TAG="latest"

FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:$DRIVER_TAG

# ARG VERSION=${VERSION:-1.9.3}

# hadolint ignore=DL3008,DL3013
RUN apt-get update && \
    apt-get install -y --no-install-recommends cmake && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    pip3 install --no-cache-dir \
        numpy==1.21.2 \
        pandas==1.3.5\
        scanpy==1.9.3
