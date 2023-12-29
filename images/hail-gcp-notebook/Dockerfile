# Ignore linting errors for pinning versions.
# hadolint ignore=DL3007
FROM gcr.io/deeplearning-platform-release/base-cpu:latest

ARG VERSION=${VERSION:-0.2.102}

COPY requirements.txt requirements.txt

# hadolint ignore=DL3008,DL3013,DL4006
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        g++ \
        liblapack3 \
        libopenblas-base \
        openjdk-8-jre-headless \
        python3-pip && \
    rm -rf /var/lib/apt/lists/* && \
    conda create -n python310 python=3.10 ipykernel && \
    conda run -n python310 pip3 install --no-cache-dir --no-deps -r requirements.txt && \
    conda run -n python310 pip3 install --no-cache-dir hail==${VERSION} && \
    curl -sSL https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-linux-x86_64.tar.bz2 | tar xjf - && \
    mv phantomjs-2.1.1-linux-x86_64/bin/phantomjs /usr/local/bin && \
    rm -r phantomjs-2.1.1-linux-x86_64 && \
    curl -o "$(conda run -n python310 find_spark_home.py)/jars/gcs-connector-hadoop2-2.0.1.jar" https://storage.googleapis.com/hadoop-lib/gcs/gcs-connector-hadoop2-2.0.1.jar
