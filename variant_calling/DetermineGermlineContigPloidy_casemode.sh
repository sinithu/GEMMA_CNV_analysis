#!/usr/bin/env bash

# https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-rare-germline-copy-number-variants
# Calling ploidies of intervals using sample read counts and ploidy model from cohort mode
# All read count files from CollectReadCounts need to be added with "-I" flags

docker run --name=determine_contig_ploidy_case -dt broadinstitute/gatk
docker exec determine_contig_ploidy_case gatk DetermineGermlineContigPloidy --model /path/contig_ploidy/ploidy-model --output-prefix ploidy-case -O /path/contig_ploidy_casemode --verbosity DEBUG -I /path/counts/sample1.hdf5 -I /path/counts/sample2.hdf5 -I ... -I /path/counts/sampleN.hdf5 2> /path/DetermineGermlineContigPloidy_casemode.log

docker stop determine_contig_ploidy_case
docker rm determine_contig_ploidy_case