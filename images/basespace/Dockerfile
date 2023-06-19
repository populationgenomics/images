FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-base:1.4

ARG VERSION=${VERSION:-1.5.3}

RUN apt-get update && \
    mkdir $HOME/bin && \
    wget "https://launch.basespace.illumina.com/CLI/${VERSION}/amd64-linux/bs" -O $HOME/bin/bs && \
    chmod u+x $HOME/bin/bs && \
    pip3 install cpg_utils && \
    pip3 install cpg_workflows && \
    mkdir $HOME/.basespace && \
    gcloud --project udn-aus secrets versions access latest --secret=basespace-config > $HOME/.basespace/basespace-config.cfg && \
    eval $($HOME/bin/bs config load basespace-config)
