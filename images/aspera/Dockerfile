FROM martinlaurent/ascli:4.18.0

USER root

# Google Cloud SDK: use the script-based installation, as the Debian package is outdated.
RUN apt update && apt install -y \
    apt-transport-https \
    bzip2 \
    ca-certificates \
    gnupg \
    skopeo \
    wget \
    zip && \
    curl https://sdk.cloud.google.com > install.sh && \
    bash install.sh --disable-prompts --install-dir=/opt && \
    rm -r install.sh /var/lib/apt/lists/* /var/cache/apt/*

USER cliuser
ENV PATH="/opt/google-cloud-sdk/bin:${PATH}"
