FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-base:1.4

SHELL ["/bin/bash", "-eo", "pipefail", "-c"]

RUN apt-get update && \
    mkdir $HOME/bin && \
    wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O $HOME/bin/bs && \
    chmod u+x $HOME/bin/bs
    