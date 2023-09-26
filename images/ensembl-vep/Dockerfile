ARG VERSION=${VERSION:-release_110.1}

FROM ensemblorg/ensembl-vep:${VERSION}

USER root

# This Dockerfile is a departure from the style used in images/vep
# Instead of building up a fresh VEP install using Miniconda, we are
# using the official VEP image and adding plugins and data files to it.
# The AlphaMissense plugin at time of writing (26.09.2023) is available
# in the vep-plugin repository, but not installed within the latest image
# the PLI, LoF, & UTR annotator data files are not available by default,
# but only add minimally to the image size.
# LOFTEE is installed by default, but the standard version doesn't work
# for GRCh38, so we overwrite the standard files
# GCSfuse and the Google Cloud SDK are not installed, as the workflow
# using this image is intended to mount data files from GCS using Hail
# Batch, so that functionality is not required inside the image.

# AlphaMissense
# print AlphaMissense scores and predictions
# only report results for the transcripts in the AlphaMissense prediction
# ./vep -i variations.vcf --plugin AlphaMissense,file=/full/path/to/file.tsv.gz,transcript_match=1
# file path = gs://dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz
ENV AlphaMissensePlugin=${VEP_DIR_PLUGINS}/AlphaMissense.pm

# pLI https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/110/pLI_values.txt
# ./vep -i variants.vcf --plugin pLI,values_file.txt
ENV PLI_SCORES=/data/pli_scores.txt

# Usage: vep -i test.vcf --tab -plugin UTRannotator,/path/to/uORF_starts_ends_GRCh38_PUBLIC.txt -o test.output
# Does not (currently) support writing  JSON output
ENV UTR38=/data/utr_annotator_38.txt

# LoFTool
# https://github.com/Ensembl/VEP_plugins/blob/release/110/LoFtool.pm
# https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/110/LoFtool_scores.txt
# ./vep -i variants.vcf --plugin LoFtool,scores_file.txt
ENV LOFTOOL_SCORES=/data/loftool_110_scores.txt

RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        git \
    && curl https://raw.githubusercontent.com/Ensembl/UTRannotator/master/uORF_5UTR_GRCh38_PUBLIC.txt > ${UTR38} \
    && curl https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/110/pLI_values.txt > ${PLI_SCORES} \
    && curl https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/110/AlphaMissense.pm > ${VEP_DIR_PLUGINS}/AlphaMissense.pm \
    && curl https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/110/LoFtool_scores.txt > ${LOFTOOL_SCORES} \
    && git clone -b grch38 https://github.com/konradjk/loftee.git \
    && mv loftee/* ${VEP_DIR_PLUGINS} \
    && rm -rf loftee \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*
