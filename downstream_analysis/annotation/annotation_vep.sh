#!/usr/bin/env bash

# Annotating CNVs using Ensembl VEP v112.0
# Download cache: https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz (https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html)
docker run --name="VEP" -dt "ensemblorg/ensembl-vep"
docker exec VEP vep -i /path/merged/samples_truvari.vcf.gz --format vcf --vcf -o /path/merged/samples_truvari_vep.vcf --fasta /path/reference.fna --offline --cache --symbol --biotype --overlaps --fork 4 --dir_cache /path/ensembl-vep/data
bgzip /path/merged/samples_truvari_vep.vcf
