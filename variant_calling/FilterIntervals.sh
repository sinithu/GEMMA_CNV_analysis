#!/usr/bin/env bash

# https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-rare-germline-copy-number-variants
# Filtering preprocessed intervals based on annotated intervals and read counts
# All read count files from CollectReadCounts need to be added with "-I" flags

docker run --name=filter_intervals -dt broadinstitute/gatk
docker exec gatk FilterIntervals -L /path/preprocessed_intervals.interval_list --annotated-intervals /path/annotated_intervals.tsv -imr OVERLAPPING_ONLY -I /path/counts/sample1.hdf5 -I /path/counts/sample2.hdf5 -I ... -I /path/counts/sampleN.hdf5 -O /path/filtered_intervals.interval_list 2> /path/FilterIntervals.log

docker stop filter_intervals
docker rm filter_intervals