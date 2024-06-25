FROM conda/miniconda3
ARG VERSION=${VERSION:-1.1.5}

RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y jasminesv==${VERSION}
