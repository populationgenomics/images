FROM australia-southeast1-docker.pkg.dev/cpg-common/images/nf_base:0.1.0-1

ENV VERSION=0.0.1

# install the latest cpg-utils code for testing the dataproc wrapper and analysis registration \
RUN git clone -b add_dataproc_handler https://github.com/populationgenomics/cpg-utils.git && \
    cd cpg-utils && \
    pip install .
