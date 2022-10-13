# Add multipy to driver-base
ARG DRIVER_TAG="latest"

FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:$DRIVER_TAG

# hadolint ignore=DL3013
RUN pip install --no-cache-dir git+https://github.com/puolival/multipy.git
