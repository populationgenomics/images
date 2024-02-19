ARG DRIVER_TAG="latest"
FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:$DRIVER_TAG
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
ARG VERSION=${VERSION:-c43c2b9c9137373c0186d3307b8089f901633128}


RUN apt-get update && \
    apt-get install -y --no-install-recommends cmake && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    # Clone the repository at a specific commit
    pip install --upgrade pip && \
    pip3 install --no-cache-dir \
        git+https://github.com/gymrek-lab/TRTools.git@$VERSION



