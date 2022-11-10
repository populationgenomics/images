# Add multipy to driver-base
ARG DRIVER_TAG="latest"

FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:$DRIVER_TAG

ARG VERSION=${VERSION:-0.16}

# hadolint ignore=DL3013
RUN pip install multipy==${VERSION}
