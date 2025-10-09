FROM us.gcr.io/broad-gatk/gatk:4.6.2.0

ARG VERSION=13.3.0

ENV EXOMISER_STEM="exomiser-cli-${VERSION}"
ENV EXOMISER_ARCHIVE="${EXOMISER_STEM}-distribution.zip"
ENV EXOMISER_URL="https://github.com/exomiser/Exomiser/releases/download/${VERSION}/${EXOMISER_ARCHIVE}"

RUN mkdir /exomiser && \
    cd /exomiser && \
    wget "${EXOMISER_URL}" && \
    unzip "${EXOMISER_ARCHIVE}" && \
    rm "${EXOMISER_ARCHIVE}" && \
    echo '#!/bin/bash' >> /usr/local/bin/exomiser && \
    echo cd /exomiser/${EXOMISER_STEM} \&\& exec java -jar $(echo /exomiser/${EXOMISER_STEM}/${EXOMISER_STEM}.jar) \"\$@\" >> /usr/local/bin/exomiser && \
    chmod 755 /usr/local/bin/exomiser

WORKDIR /exomiser/${EXOMISER_STEM}

# requires at least Java 17
# download & unzip the distribution - this will create a directory called 'exomiser-cli-${VERSION}'
RUN apt-get update && \
    apt-get install --no-install-recommends -y \
		python3 \
		python3-pip \
		tree && \
    apt-get clean; \
    pip3 install setuptools && \
    pip3 install PyYAML

COPY config_shuffle.py template.yaml ./
