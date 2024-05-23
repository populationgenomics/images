ARG DRIVER_TAG="latest"

FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:$DRIVER_TAG

# hadolint ignore=DL3008,DL3013
RUN apt-get update && \
    apt-get install -y --no-install-recommends cmake && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    pip3 install --no-cache-dir \
        pandas==2.2.1 \
        openpyxl==3.1.2