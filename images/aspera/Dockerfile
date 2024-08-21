FROM google/cloud-sdk:486.0.0

# Aspera doesn't install for the root user.
RUN useradd -m aspera
USER aspera
WORKDIR /home/aspera

RUN curl https://d3gcli72yxqn2z.cloudfront.net/downloads/connect/latest/bin/ibm-aspera-connect_4.2.12.780_linux_x86_64.tar.gz | tar xz && \
    ./ibm-aspera-connect_4.2.12.780_linux_x86_64.sh && \
    rm ibm-aspera-connect_4.2.12.780_linux_x86_64.sh

# Allow Hail Batch to create files in /io.
USER root
