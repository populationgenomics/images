from debian:10.13

ARG VERSION=${VERSION}

RUN apt-get update && apt-get install --yes \
    build-essential \
    gcc-multilib \
    gengetopt \
    libgcc1 \
    libc6 \
    liblapack-dev \
    libpthread-stubs0-dev \
    libstdc++6 \
    procps \
    unzip \
    wget \ 
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/local/

RUN wget https://bitbucket.org/statgen/xhmm/get/master.zip && \
    unzip master.zip && \
    rm master.zip

WORKDIR /usr/local/statgen-xhmm-998f7c405974

RUN make && \
    ln -s $(realpath build/execs/xhmm)  /usr/local/bin/xhmm

CMD []
