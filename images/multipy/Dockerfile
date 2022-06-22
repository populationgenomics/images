# Add multipy to driver-base
ARG DRIVER_SHA="36c6d4548ef347f14fd34a5b58908057effcde82-hail-83056327f288917537531475ba475287b413db1c"

FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:$DRIVER_SHA

# hadolint ignore=DL3013
RUN pip install --no-cache-dir git+https://github.com/puolival/multipy.git
