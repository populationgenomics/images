FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:a15ed15c0ca0b7304e0ceca1b873d734720c5d7a-hail-abaf45c99b32250dfc55a1b26a6a0635760e0a99

# ARG VERSION=${VERSION:-0.2.124}

# hadolint ignore=DL3008,DL3013
RUN apt-get update && \
    apt-get install -y --no-install-recommends cmake && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    pip3 install --no-cache-dir \
        hail==0.2.124