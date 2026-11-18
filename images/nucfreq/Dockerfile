FROM python:3.11-slim-bookworm

# NucFreq has no releases or tags, so the version is the git commit of mrvollger/NucFreq that
# NucPlot.py is taken from (2025-10-06, the PDF rasterising merge).
ARG VERSION=88ac15a819b2fdbd81806d1a5123e1d8ba03f674

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# procps provides ps, which Nextflow needs to record task CPU and memory use.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        procps \
        wget && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

# NucPlot.py predates numpy 2 and pandas 2 (it uses DataFrame.append); pysam reads the BAM, the rest draw.
RUN pip install --no-cache-dir \
        "matplotlib==3.10.8" \
        "numpy==1.26.4" \
        "pandas==1.5.3" \
        "pysam==0.23.3" \
        "seaborn==0.13.2" && \
    wget -qO /usr/local/bin/NucPlot.py "https://raw.githubusercontent.com/mrvollger/NucFreq/${VERSION}/NucPlot.py" && \
    chmod +x /usr/local/bin/NucPlot.py

# The script parses its arguments and imports every library before touching a BAM.
RUN NucPlot.py --help > /dev/null
