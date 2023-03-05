FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-base:1.4

# Override this with a specific SHA in images.toml, as haplocheck doesn't have releases.
ARG VERSION=HEAD

    git clone --depth 1 https://github.com/leklab/haplocheckCLI.git && \
    cd haplocheckCLI && \
    git checkout $VERSION && \
    cp haplocheckCLI.jar ../ && \
    cd .. && \
    rm -rf haplocheckCLI
