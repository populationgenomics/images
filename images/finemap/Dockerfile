FROM debian:bullseye-slim

WORKDIR /tmp

# Install dev libraries
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    curl \
    libgomp1 \
    && apt-get clean

RUN curl -L http://www.christianbenner.com/finemap_v1.4.2_x86_64.tgz -o finemap.tgz \
    && tar zxvf finemap.tgz && mv finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64 /usr/local/bin/finemap \
    && chmod +x /usr/local/bin/finemap && rm -rf /tmp/*

CMD ["bash"]
