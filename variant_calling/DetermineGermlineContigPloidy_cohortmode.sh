#!/usr/bin/env bash

# https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-rare-germline-copy-number-variants
# Calling ploidies of intervals using sample read counts and contig ploidy priors table (with prior ploidy probabilities of chromosomes)
# All read count files from CollectReadCounts need to be added with "-I" flags

docker run --name=determine_contig_ploidy_cohort -dt broadinstitute/gatk
docker exec determine_contig_ploidy_cohort gatk DetermineGermlineContigPloidy -L /path/filtered_intervals.interval_list -imr OVERLAPPING_ONLY --contig-ploidy-priors /path/ploidy_priors.tsv -O /path/contig_ploidy --output-prefix ploidy --verbosity DEBUG -I /path/counts/sample1.hdf5 -I /path/counts/sample2.hdf5 -I ... -I /path/counts/sampleN.hdf5 2> /path/DetermineGermlineContigPloidy_cohortmode.log

docker stop determine_contig_ploidy_cohort
docker rm determine_contig_ploidy_cohort