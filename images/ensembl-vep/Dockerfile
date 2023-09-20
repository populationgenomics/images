ARG VERSION=${VERSION:-release_110.1}

FROM ensemblorg/ensembl-vep:${VERSION}

USER root

## Copying in a couple of small data files, and a plugin or two

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
    #     apt-transport-https \
    #
        ca-certificates \
        git \
    #     gnupg \
    #     lsb-release \

    ## I don't think we need google-cloud-sdk inside the image..
    # && echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" \
    #| tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
    #| apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - \
    #&& apt-get update \
    #&& apt-get install google-cloud-cli -y \
    && curl https://raw.githubusercontent.com/Ensembl/UTRannotator/master/uORF_5UTR_GRCh38_PUBLIC.txt > ${UTR38} \
    && curl https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/110/pLI_values.txt > ${PLI_SCORES} \
    && curl https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/110/AlphaMissense.pm > ${VEP_DIR_PLUGINS}/AlphaMissense.pm \
    && curl https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/110/LoFtool_scores.txt > ${LOFTOOL_SCORES} \
    ## LOFTEE files - main branch of LOFTEE is GRCh37 only
    && git clone -b grch38 https://github.com/konradjk/loftee.git \
    && cp loftee/*.p* ${VEP_DIR_PLUGINS} \
    && rm -rf loftee \

    ## Install gcsfuse to mount VEP cache files
    ## Do we need to have gcsfuse usable from inside the container? It's mounted on by Hail Batch
    #&& GCSFUSE_REPO=gcsfuse-$(lsb_release -c -s) \
    #&& echo "deb http://packages.cloud.google.com/apt $GCSFUSE_REPO main" | tee /etc/apt/sources.list.d/gcsfuse.list \
    #&& curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - \
    #&& apt-get update \
    #&& apt-get install -y gcsfuse \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*
