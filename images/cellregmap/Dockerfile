ARG DRIVER_TAG="latest"

FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:$DRIVER_TAG

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# hadolint ignore=DL3008,DL3013
RUN apt-get update && \
    apt-get install -y --no-install-recommends cmake && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    # The chi2comb-py Python package depends on the chi2comb C++ library.
    curl -fsSL https://raw.githubusercontent.com/limix/chi2comb/master/install | bash && \
    pip3 install --no-cache-dir \
        cellregmap==0.0.3 \
        chi2comb==0.1.0 \
        click==8.0.4 \
        limix==3.0.4 \
        numpy==1.21.2 \
        pandas==1.3.4 \
        pandas_plink==2.2.9 \
        scanpy==1.8.2 \
        scipy==1.7.3 \
        xarray==0.20.2
