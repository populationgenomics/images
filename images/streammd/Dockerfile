ARG ST_VERSION=1.21-1

FROM australia-southeast1-docker.pkg.dev/cpg-common/images/samtools:${ST_VERSION}

# the streammd version
ARG VERSION=4.3.0

# Install the runtime shared libraries that the compiled htslib may depend on
RUN apt-get update && apt-get install -y \
		cmake \
		g++ \
		gcc \
		git \
		make \
		wget \
    && rm -rf /var/lib/apt/lists/* && \
    wget https://github.com/delocalizer/streammd/releases/download/${VERSION}/streammd-${VERSION}.tar.gz && \
		tar -xzvf streammd-${VERSION}.tar.gz && \
		cd streammd-${VERSION} && \
		./configure && \
		make && \
		make check && \
		make install && \
		cd .. && \
		rm -r streammd-${VERSION} && \
    streammd --help
