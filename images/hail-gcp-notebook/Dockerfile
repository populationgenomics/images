FROM gcr.io/deeplearning-platform-release/base-cpu

RUN apt-get update && \
    apt-get install -y --no-install-recommends openjdk-8-jre-headless g++ python3-pip libopenblas-base liblapack3 && \
    rm -rf /var/lib/apt/lists/* && \
    pip3 install --no-cache-dir "cloudpathlib[all]" hail phantomjs pyarrow && \
    cd $(find_spark_home.py)/jars && \
    curl -O https://storage.googleapis.com/hadoop-lib/gcs/gcs-connector-hadoop2-2.0.1.jar
