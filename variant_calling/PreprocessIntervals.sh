#!/usr/bin/env bash

# https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-rare-germline-copy-number-variants
# Didiving reference genome (GRCh38/hg38) to intervals of 1000 bp

# Excluding ENCODE blacklisted genomic regions, download: https://github.com/Boyle-Lab/Blacklist/tree/master/lists/Blacklist_v1 (publication: https://doi.org/10.1038/s41598-019-45839-z)

docker run --name=preprocess -dt broadinstitute/gatk
docker exec preprocess gatk PreprocessIntervals -XL /path/hg38-blacklist.bed -R /path/reference.fna --padding 0 -imr OVERLAPPING_ONLY -O /path/preprocessed_intervals.interval_list 2> /path/PreprocessIntervals.log

docker stop preprocess
docker rm preprocess
