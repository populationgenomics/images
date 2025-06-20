FROM australia-southeast1-docker.pkg.dev/cpg-common/images/dragmap:1.3.0-tokenizer-next-fix

ARG VERSION="2.7.0"

# source
# https://sapac.support.illumina.com/downloads/dragen_decompression_software-linux.html

# installation
# https://support-docs.illumina.com/SW/Decompression_v2_6/Content/SW/ORA/Decompression_Installation.htm

RUN wget -q https://s3.amazonaws.com/webdata.illumina.com/downloads/software/dragen-decompression/orad.2.7.0.linux.tar.gz && \
    tar -xzf orad.2.7.0.linux.tar.gz && \
    mv orad.2.7.0.linux/orad /usr/local/bin/orad && \
    rm -rf orad.2.7.0.linux.tar.gz orad.2.7.0.linux
