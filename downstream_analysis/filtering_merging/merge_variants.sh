#!/usr/bin/env bash
# Merging CNVs that are close to each other
truvari collapse -i /path/merged/samples.vcf.gz -o /path/merged/samples_truvari.vcf --pctseq 0 --pctsize=0.5 --refdist=1000
bgzip /path/merged/samples_truvari.vcf
tabix /path/merged/samples_truvari.vcf.gz