FROM gcr.io/deeplearning-platform-release/base-cpu:m93

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        g++ \
        liblapack3 \
        libopenblas-base \
        openjdk-8-jre-headless \
        python3-pip && \
    rm -rf /var/lib/apt/lists/* && \
    pip3 install --no-cache-dir \
        "cloudpathlib[all]==0.9.0" \
        cpg-utils==4.4.3 \
        hail==0.2.96 \
        phantomjs==1.4.1 \
        pyarrow==8.0.0 \
        sample-metadata==4.11.0 && \
    curl -o "$(find_spark_home.py)/jars/gcs-connector-hadoop2-2.0.1.jar" https://storage.googleapis.com/hadoop-lib/gcs/gcs-connector-hadoop2-2.0.1.jar
