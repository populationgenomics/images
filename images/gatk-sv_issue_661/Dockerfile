FROM australia-southeast1-docker.pkg.dev/cpg-common/images/sv/sv-pipeline:2023-09-13-v0.28.3-beta-af8362e3

ARG VERSION=${VERSION:-n_a}

COPY compute_AFs.py /opt/sv-pipeline/05_annotation/scripts/compute_AFs.py

RUN chmod +x /opt/sv-pipeline/05_annotation/scripts/compute_AFs.py
