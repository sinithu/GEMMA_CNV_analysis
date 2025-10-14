#!/usr/bin/env bash

# https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-rare-germline-copy-number-variants
# Annotating intervals with GC and mappability content
# Mappability track: Umap single-read k100, download https://bismap.hoffmanlab.org (publication: https://doi.org/10.1093/nar/gky677)

docker run --name=annotate_intervals -dt broadinstitute/gatk
docker exec annotate_intervals gatk AnnotateIntervals -R /path/reference.fna -L /path/preprocessed_intervals.interval_list --mappability-track /path/merged_k100.umap.bed.gz -imr OVERLAPPING_ONLY -O /path/annotated_intervals.tsv 2> /path/AnnotateIntervals.log

docker stop annotate_intervals
docker rm annotate_intervals