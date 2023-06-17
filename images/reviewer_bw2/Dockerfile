FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

ARG VERSION=${VERSION:-0.2.7}

# Install Git to clone the repository
RUN apt-get update && apt-get install -y git cmake g++ zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev

# Clone the repository
RUN git clone https://github.com/bw2/REViewer

# Set the working directory to the cloned repository
WORKDIR /REViewer

# Create a build directory and navigate to it
RUN mkdir build
WORKDIR /REViewer/build

# Run CMake to configure the project
RUN cmake ..

# Build the executable
RUN make

# Debugging steps
RUN ls -l /REViewer/build    # Print the contents of the build directory

# Set the entry point for the Docker container
ENTRYPOINT ["/REViewer/build"]

