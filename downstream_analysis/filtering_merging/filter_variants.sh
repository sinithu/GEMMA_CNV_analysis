#!/usr/bin/env bash

#Adding symbolic links of VCFs from case and cohort mode to the same directory
mkdir /path/postprocess
ln -s /path/postprocess_case/* /path/postprocess
ln -s /path/postprocess_cohort/* /path/postprocess

#Filtering variants in segments-VCFs based on qualities
cd /path/postprocess
files=$(ls *segments.vcf.gz)
mkdir /path/filtered_VCFs

for file in $files
do
	bcftools view -h $file > /path/filtered_VCFs/${file%.vcf.gz}_filtered.vcf
	bcftools view -H -i 'QA>=20 | QS>=20 | QSE>=20 | QSS>=20' $file >> /path/filtered_VCFs/${file%.vcf.gz}_filtered.vcf
	bgzip /path/filtered_VCFs/${file%.vcf.gz}_filtered.vcf
	tabix /path/filtered_VCFs/${file%.vcf.gz}_filtered.vcf.gz
done
