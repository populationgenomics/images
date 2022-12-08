# Add multipy to driver-base
ARG DRIVER_TAG="latest"

FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:$DRIVER_TAG

# Install a Python3-compatible version that hasn't been released on PyPI yet.
RUN pip install git+https://github.com/puolival/multipy.git@ef9d2e2ffd6a7e079e29268dfbd59e630dd346b4

