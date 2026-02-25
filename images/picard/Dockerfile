# --- STAGE 1: Build Picard from Source ---
FROM eclipse-temurin:17-jdk AS builder

ENV MAMBA_ROOT_PREFIX /root/micromamba
ENV PATH $MAMBA_ROOT_PREFIX/bin:$PATH
ARG VERSION=3.4.0

RUN apt-get update && \
    apt-get --no-install-recommends install -y git && \
    rm -rf /var/lib/apt/lists/*

# Clone Picard repository
RUN git clone --branch 3.4.0 --depth 1 https://github.com/broadinstitute/picard.git /usr/picard
WORKDIR /usr/picard

# Build Picard with htsjdk 4.3.0 for CRAM v3.1 support
RUN ./gradlew clean printVersion shadowJar -Dhtsjdk.version=4.3.0

# --- STAGE 2: Minimal Runtime Environment ---
FROM debian:bookworm-slim
ENV MAMBA_ROOT_PREFIX /root/micromamba
ENV PATH $MAMBA_ROOT_PREFIX/bin:$PATH
ARG VERSION=3.4.0

# Install minimal system requirements for Micromamba
RUN apt-get update && apt-get install -y wget bzip2 ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Install Micromamba
RUN wget -qO- https://micro.mamba.pm/api/micromamba/linux-64/latest | \
    tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba && \
    chmod +x /usr/local/bin/micromamba

# Install samtools, gcloud-sdk, and the Java Runtime (JRE)
# We use 'openjdk' here because the JAR needs a Java environment to run.
RUN mkdir -p ${MAMBA_ROOT_PREFIX} && \
    micromamba install -y --prefix ${MAMBA_ROOT_PREFIX} \
    -c bioconda -c conda-forge \
    samtools \
    google-cloud-sdk \
    openjdk=17 && \
    micromamba clean --all --yes

# Copy the JAR from the builder stage
RUN mkdir /usr/picard
COPY --from=builder /usr/picard/build/libs/picard.jar /usr/picard/picard.jar

# Setup working directory
WORKDIR /usr/working
