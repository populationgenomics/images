FROM debian:bookworm-slim AS installer

ARG VERSION=2.36.29

RUN apt-get update && apt-get install -y \
        curl \
        unzip && \
    curl -s "https://awscli.amazonaws.com/awscli-exe-linux-$(uname -m)-${VERSION}.zip" -o awscliv2.zip && \
    unzip -q awscliv2.zip && \
    ./aws/install --bin-dir /usr/local/bin --install-dir /usr/local/aws-cli

FROM debian:bookworm-slim

ARG VERSION

LABEL \
  version="v${VERSION}" \
  description="AWS CLI v2"

# groff-base and less are required by "aws help"
RUN apt-get update && apt-get install -y \
        ca-certificates \
        groff-base \
        less && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

COPY --from=installer /usr/local/aws-cli /usr/local/aws-cli

RUN ln -s /usr/local/aws-cli/v2/current/bin/aws /usr/local/bin/aws && \
    ln -s /usr/local/aws-cli/v2/current/bin/aws_completer /usr/local/bin/aws_completer
