FROM ubuntu:18.04

ARG VERSION=${VERSION:-v0.3.15}
ARG HTS_VERSION=1.15.1

# installation from https://github.com/Illumina/hap.py/blob/master/Dockerfile
ENV DEBIAN_FRONTEND=noninteractive

# hadolint ignore=DL3008
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        autoconf \
        build-essential \
        bzip2 \
        cmake \
        cython \
        git \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        openjdk-8-jdk \
        pkg-config \
        python \
        python2.7 \
        python2.7-dev \
        python-setuptools \
        python-pip \
        python-psutil \
        python-numpy \
        python-pandas \
        python-distribute \
        python-pysam \
        python-scipy \
        software-properties-common \
        wget \
        zlib1g-dev && \
    apt-get clean -y  && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir bx-python==0.8.9

# get + install ant
WORKDIR /opt
RUN wget -q http://archive.apache.org/dist/ant/binaries/apache-ant-1.9.7-bin.tar.gz && \
    tar xzf apache-ant-1.9.7-bin.tar.gz && \
    rm apache-ant-1.9.7-bin.tar.gz
ENV PATH="$PATH:/opt/apache-ant-1.9.7/bin"

# clone in, install, delete - then update bcftools/samtools
RUN git clone --depth 1 --single-branch --branch ${VERSION} https://github.com/Illumina/hap.py /opt/hap.py-source && \
	python /opt/hap.py-source/install.py /opt/hap.py --with-rtgtools --no-tests && \
	rm -rf /opt/hap.py-source && \
	wget https://github.com/samtools/bcftools/releases/download/${HTS_VERSION}/bcftools-${HTS_VERSION}.tar.bz2 && \
	tar -xf bcftools-${HTS_VERSION}.tar.bz2 && \
	cd bcftools-${HTS_VERSION} && \
	bash configure --prefix=/opt/hap.py && \
	make && \
	make install && \
	cd .. && \
	rm -rf bcftools-${HTS_VERSION}.tar.bz2 bcftools-${HTS_VERSION} && \
	wget https://github.com/samtools/samtools/releases/download/${HTS_VERSION}/samtools-${HTS_VERSION}.tar.bz2 && \
	tar -xf samtools-${HTS_VERSION}.tar.bz2 && \
	cd samtools-${HTS_VERSION} && \
	bash configure --prefix=/opt/hap.py && \
	make && \
	make install && \
	cd .. && \
	rm -rf samtools-${HTS_VERSION}.tar.bz2 samtools-${HTS_VERSION}

ENV PATH="/opt/hap.py/bin:$PATH"
ENV RTG_JAR="/opt/hap.py/libexec/rtg-tools-install/RTG.jar"
# run the RTG JAR file with $ java -jar "RTG_JAR"
