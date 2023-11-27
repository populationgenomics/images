FROM r-base:4.3.1

ARG VERSION=${VERSION}

LABEL software="ExomeDepth"
LABEL description="R package ExomeDepth is used to perform CNV calling on a library of BAM files"
LABEL documentation="https://CRAN.R-project.org/package=ExomeDepth"

RUN apt-get -y update && \
	apt-get --yes --force-yes install \
        libcurl4-openssl-dev \ 
        procps 

RUN R -e "install.packages('BiocManager'); \
    BiocManager::install(); \
    BiocManager::install(c( \
        'Biostrings', \
        'GenomicAlignments', \
        'GenomicRanges', \
        'IRanges', \
        'Rsamtools' \
    )); \
    install.packages(c( \
        'ExomeDepth', \
        'optparse' \
    ))"

RUN mkdir /data && \
	mkdir /data/rscripts/

COPY . /data/rscripts/

WORKDIR /data/rscripts/

