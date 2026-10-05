# Use a slim Python image for efficiency
FROM python:3.11-slim

ARG VERSION=7.14.0

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*

# Install metamist and any other dependencies
# (Nextflow's bin/ script is automatically added to the PATH by Nextflow,
#  so we don't strictly need to COPY it into the container if we mount it,
#  but for a truly self-contained container on GCP AR, we can include it.)
RUN pip install --no-cache-dir "metamist>=7.14.0"

# (Optional) If you want the script inside the image:
# COPY bin/fetch_metamist.py /usr/local/bin/fetch_metamist.py
# RUN chmod +x /usr/local/bin/fetch_metamist.py

# Standard entrypoint (though Nextflow overrides this)
CMD ["python3"]
