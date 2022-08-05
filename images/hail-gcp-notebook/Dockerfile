# Ignore linting errors for pinning versions.
# hadolint ignore=DL3007
FROM gcr.io/deeplearning-platform-release/base-cpu:latest

# hadolint ignore=DL3008,DL3013
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        g++ \
        liblapack3 \
        libopenblas-base \
        openjdk-8-jre-headless \
        python3-pip && \
    rm -rf /var/lib/apt/lists/* && \
    pip3 install --no-cache-dir \
        bokeh \
        "cloudpathlib[all]" \
        cpg-utils \
        hail \
        phantomjs \
        pyarrow \
        sample-metadata \
        selenium && \
    curl -o "$(find_spark_home.py)/jars/gcs-connector-hadoop2-2.0.1.jar" https://storage.googleapis.com/hadoop-lib/gcs/gcs-connector-hadoop2-2.0.1.jar
