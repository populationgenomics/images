ARG DRIVER_TAG="e2aa072191014ba1db0a02eb5e868a4daf327436-hail-d85618aaf4d2bf667aa72c5c6a142e9bafaae9c8"

FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:$DRIVER_TAG

ARG VERSION=${VERSION:-3.5.1}
RUN pip install --no-cache-dir paramiko==${VERSION}
