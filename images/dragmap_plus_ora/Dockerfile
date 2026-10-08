FROM australia-southeast1-docker.pkg.dev/cpg-common/images/dragmap:post-bullseye-EOL-fixes-1

ARG VERSION=2.7.0

# source
# https://sapac.support.illumina.com/downloads/dragen_decompression_software-linux.html

# installation
# https://support-docs.illumina.com/SW/Decompression_v2_6/Content/SW/ORA/Decompression_Installation.htm

RUN wget -q https://s3.amazonaws.com/webdata.illumina.com/downloads/software/dragen-decompression/orad.${VERSION}.linux.tar.gz && \
    tar -xzf orad.${VERSION}.linux.tar.gz && \
    mv orad.${VERSION}.linux/orad /usr/local/bin/orad && \
    rm -rf orad.${VERSION}.linux.tar.gz orad.${VERSION}.linux
