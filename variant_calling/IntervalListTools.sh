#!/usr/bin/env bash

# https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-rare-germline-copy-number-variants
# Scattering interval list into multiple lists

docker run --name=interval_list_tools -dt broadinstitute/gatk
docker exec interval_list_tools gatk IntervalListTools --INPUT /path/filtered_intervals.interval_list --SUBDIVISION_MODE INTERVAL_COUNT --SCATTER_CONTENT 5000 --OUTPUT /path/scatter 2> /path/IntervalListTools.log

docker stop interval_list_tools
docker rm interval_list_tools