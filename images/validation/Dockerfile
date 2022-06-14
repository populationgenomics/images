FROM pkrusche/hap.py:v0.3.9

# hadolint ignore=DL3008
RUN apt-get update && \
    apt-get install -y unzip && \
    apt-get clean -y

# Install VcfEval
RUN wget https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.12.1/rtg-tools-3.12.1-linux-x64.zip && \
    unzip rtg-tools-3.12.1-linux-x64.zip && \
    rm rtg-tools-3.12.1-linux-x64.zip && \
    mv rtg-tools-3.12.1 /vcfeval

ENV PATH="/opt/hap.py/bin:$PATH" VCFEVAL_JAR="/vcfeval/RTG.jar"
