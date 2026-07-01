# Use a slim Python image for efficiency
FROM python:3.12-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*

# Install metamist and any other dependencies
ARG VERSION=${VERSION:-2.4.1}
RUN pip install --no-cache-dir "rbceq2==${VERSION}"

# Standard entrypoint (though Nextflow overrides this)
CMD ["python3"]
