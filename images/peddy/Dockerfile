# This image adds cpg-hail, SM and peddy to the base image.
FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-base:1.0

ENV HAIL_QUERY_BACKEND service

RUN apt update && apt install -y g++ libopenblas-base liblapack3 openjdk-11-jre-headless python3-pip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    pip3 install \
        analysis-runner \
        bokeh \
        cpg-hail \
        cpg-utils \
        peddy \
        sample-metadata 