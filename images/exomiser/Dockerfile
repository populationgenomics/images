FROM us.gcr.io/broad-gatk/gatk:4.4.0.0
ARG VERSION=${VERSION:-13.3.0}
ARG DATA_VERSION=${DATA_VERSION:-2302}
ENV EXOMISER_STEM="exomiser-cli-${VERSION}"
ENV EXOMISER_ARCHIVE="${EXOMISER_STEM}-distribution.zip"
ENV EXOMISER_URL="https://github.com/exomiser/Exomiser/releases/download/${VERSION}/${EXOMISER_ARCHIVE}"
RUN mkdir /exomiser && \
    cd /exomiser && \
    wget "${EXOMISER_URL}" && \
    unzip "${EXOMISER_ARCHIVE}" && \
    rm "${EXOMISER_ARCHIVE}" && \
    echo '#!/bin/bash' >> /usr/local/bin/exomiser && \
    echo cd /exomiser/${EXOMISER_STEM} \&\& exec java -jar /exomiser/${EXOMISER_STEM}/${EXOMISER_STEM}.jar \"\$@\" >> /usr/local/bin/exomiser && \
    chmod 755 /usr/local/bin/exomiser

WORKDIR /exomiser/${EXOMISER_STEM}

# download & unzip the distribution - this will create a directory called 'exomiser-cli-14.0.0'
RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    python3 \
    python3-pip \
    tree && \
    apt-get clean; \
    pip3 install setuptools && \
    pip3 install PyYAML

COPY config_shuffle.py template.yaml ./

# create the application.properties file, responsive to version ARGs
RUN echo "exomiser.data-directory=/exomiser/${EXOMISER_STEM}/data" > /exomiser/${EXOMISER_STEM}/application.properties && \
    echo "exomiser.hg38.data-version=${DATA_VERSION}" >> /exomiser/${EXOMISER_STEM}/application.properties && \
    echo "exomiser.phenotype.data-version=${DATA_VERSION}" >> /exomiser/${EXOMISER_STEM}/application.properties && \
    echo "exomiser.hg38.clin-var-data-version=${DATA_VERSION}" >> /exomiser/${EXOMISER_STEM}/application.properties && \
    echo "exomiser.hg38.use-clinvar-white-list=true" >> /exomiser/${EXOMISER_STEM}/application.properties && \
    echo "logging.level.com.zaxxer.hikari=ERROR" >> /exomiser/${EXOMISER_STEM}/application.properties
