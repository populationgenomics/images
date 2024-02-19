ARG DRIVER_TAG="latest"
FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:$DRIVER_TAG
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
ARG VERSION=${VERSION:-1d61fd1011dd3f6e0cf1c41cb21a32b8b8be208c}


RUN apt-get update && \
    apt-get install -y --no-install-recommends cmake && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    # Clone the repository at a specific commit
    pip3 install --upgrade pip && \
    pip3 install --no-cache-dir \
        git+https://github.com/gymrek-lab/TRTools.git@$VERSION