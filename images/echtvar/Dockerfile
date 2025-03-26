FROM debian:bookworm-slim

ARG VERSION=${VERSION:-v0.2.1}

ADD https://github.com/brentp/echtvar/releases/download/${VERSION}/echtvar /bin/echtvar

RUN chmod +x /bin/echtvar

# designed for the gnomAD v4 joint-dataset
COPY config.json /config.json

# designed for the gnomAD v4 genomes-only
COPY simple_config.json /simple_config.json

ENV ECHTVAR_CONFIG=/config.json
ENV ECHTVAR_SIMPLE_CONFIG=/simple_config.json
