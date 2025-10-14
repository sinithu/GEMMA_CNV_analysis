# GEMMA CNV analysis
### Pipeline for copy number variation calling and analysis using preprocessed GEMMA WGS data
Germline copy number variations (CNVs) were called using GATK-gCNV tools (https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-rare-germline-copy-number-variants). Cohort mode was performed using PreprocessIntervals, AnnotateIntervals, CollectReadCounts, FilterIntervals, DetermineGermlineContigPloidy, IntervalListTools, GermlineCNVCaller and PostprocessGermlineCNVCalls with input BAM files and reference fasta file. Case mode was performed with case-mode specific runs using CollectReadCounts, DetermineGermlineContigPloidy, GermlineCNVCaller and PostprocessGermlineCNVCalls. 

Produced VCF files were filtered based on variant qualities and merged to multisample VCF using bcftools view and truvari collapse. CNVs were annotated using Ensembl VEP.

![alt text](https://github.com/user-attachments/files/22881782/gemma_cnv_pipeline.pdf)

Scripts for GATK germline CNV variant calling can be found in [variant_calling](https://github.com/sinithu/GEMMA_CNV_analysis/tree/main/variant_calling).
Scripts for downstream analyses can be found in [downstream_analysis](https://github.com/sinithu/GEMMA_CNV_analysis/tree/main/downstream_analysis).
