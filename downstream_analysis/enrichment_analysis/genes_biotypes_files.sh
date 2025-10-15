#!/usr/bin/env bash

#Dividing annotated multisample VCF to cases and controls using sample lists
#Both files contain the same amount of rows (also the rows without variant in any of the samples)
bcftools view -S /path/cases.txt /path/merged/samples_truvari_vep.vcf.gz -Oz -o samples_truvari_vep_cases.vcf.gz
bcftools view -S /path/controls.txt /path/merged/samples_truvari_vep.vcf.gz -Oz -o samples_truvari_vep_controls.vcf.gz

#Creating files of gene information (VEP columns "SYMBOL" and "BIOTYPE") separately in deletions and duplications
#COUNT(CN!=".")>0 is used to include only actual variants (i.e. variants with copy number information in at least one sample)
bcftools view -i 'ALT="<DEL>" & COUNT(CN!=".")>0' samples_truvari_vep_cases.vcf.gz | bcftools +split-vep -f'%SYMBOL\t%BIOTYPE\n' > genes_biotype_DEL_cases.txt
bcftools view -i 'ALT="<DUP>" & COUNT(CN!=".")>0' samples_truvari_vep_cases.vcf.gz | bcftools +split-vep -f'%SYMBOL\t%BIOTYPE\n' > genes_biotype_DUP_cases.txt
bcftools view -i 'ALT="<DEL>" & COUNT(CN!=".")>0' samples_truvari_vep_controls.vcf.gz | bcftools +split-vep -f'%SYMBOL\t%BIOTYPE\n' > genes_biotype_DEL_controls.txt
bcftools view -i 'ALT="<DUP>" & COUNT(CN!=".")>0' samples_truvari_vep_controls.vcf.gz | bcftools +split-vep -f'%SYMBOL\t%BIOTYPE\n' > genes_biotype_DUP_controls.txt
