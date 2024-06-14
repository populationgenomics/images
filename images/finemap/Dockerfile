FROM debian:stretch-slim

WORKDIR /tmp

# Install dev libraries
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    curl \
    unzip \
    gzip \
    bzip2 \
    ca-certificates \
    build-essential \
    gfortran \
    libgfortran-6-dev \
    libgomp1 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

# Fine-mapping tools required packages
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    libgsl-dev libgsl2 libatlas3-base liblapack-dev \
    && apt-get clean

# R, Python, SoS and DSC
ENV MINICONDA_VERSION 4.7.12
ENV PATH /opt/miniconda3/bin:$PATH

RUN curl https://repo.anaconda.com/miniconda/Miniconda3-$MINICONDA_VERSION-Linux-x86_64.sh -o MCON.sh \
    && /bin/bash MCON.sh -b -p /opt/miniconda3 \
    && ln -s /opt/miniconda3/etc/profile.d/conda.sh /etc/profile.d/conda.sh

RUN curl -L http://www.christianbenner.com/finemap_v1.4.2_x86_64.tgz -o finemap.tgz \
    && tar zxvf finemap.tgz && mv finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64 /usr/local/bin/finemap \
    && chmod +x /usr/local/bin/finemap && rm -rf /tmp/*

CMD ["bash"]