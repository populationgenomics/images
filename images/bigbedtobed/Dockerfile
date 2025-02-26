FROM debian:bookworm-slim

ARG VERSION=${VERSION:-1.04.00}

ADD http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed /usr/local/bin/

RUN chmod a+x /usr/local/bin/bigBedToBed
