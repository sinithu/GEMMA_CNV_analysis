#!/usr/bin/env bash
# Merging filtered VCFs using bcftools merge
# Add the all samples to the end of the command
bcftools merge -m none -O z -o /path/merged/samples.vcf.gz /path/filtered_VCFs/sample1.vcf.gz /path/filtered_VCFs/sample2.vcf.gz ... /path/filtered_VCFs/sampleN.vcf.gz