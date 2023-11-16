FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-base:1.4

ARG VERSION=${VERSION:-1.5.3}

RUN apt-get update && \
    mkdir $HOME/bin && \
    wget "https://launch.basespace.illumina.com/CLI/${VERSION}/amd64-linux/bs" -O $HOME/bin/bs && \
    chmod u+x $HOME/bin/bs && \
    mkdir $HOME/.basespace
