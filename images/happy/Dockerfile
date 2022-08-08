FROM ubuntu:18.04

# installation from https://github.com/Illumina/hap.py/blob/master/Dockerfile
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y \
        autoconf \
        build-essential \
        bzip2 \
        cmake \
        cython \
        git \
        libbz2-dev \
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
    apt-get clean -y

RUN pip install bx-python

# get + install ant
WORKDIR /opt
RUN wget http://archive.apache.org/dist/ant/binaries/apache-ant-1.9.7-bin.tar.gz && \
    tar xzf apache-ant-1.9.7-bin.tar.gz && \
    rm apache-ant-1.9.7-bin.tar.gz
ENV PATH $PATH:/opt/apache-ant-1.9.7/bin

# clone in, install, delete
RUN git clone --depth 1 --branch v0.3.15 https://github.com/Illumina/hap.py /opt/hap.py-source && \
	python /opt/hap.py-source/install.py /opt/hap.py --with-rtgtools --no-tests && \
	rm -rf /opt/hap.py-source

ENV PATH="/opt/hap.py/bin:$PATH"
ENV RTG_JAR="/opt/hap.py/libexec/rtg-tools-install/RTG.jar"
# run the RTG JAR file with $ java -jar "$VCFEVAL_JAR"
