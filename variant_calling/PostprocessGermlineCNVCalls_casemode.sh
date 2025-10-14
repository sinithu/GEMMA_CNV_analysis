#!/usr/bin/env bash

# https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-rare-germline-copy-number-variants
# For CASE MODE
# Creating VCF files of copy number calls
# This script first outputs a run file with PostprocessGermlineCNVCalls commands for all the samples, and then executes it

#Creating a PostprocessGermlineCNVCalls file for all samples and scatters
#All calls paths from GermlineCNVCaller (case mode) and and model paths from GermlineCNVCaller (cohort mode) need to be added with "--calls-shard-path" and "--model-shard-path" flags
#For 72 samples
export run_file="/path/postprocessCNV_casemode_samples.sh"
echo "#!/usr/bin/env bash" > $run_file
mkdir /path/postprocess_case

for i in {0..71}
do
    cd "/path/case_gCNV_calls/case_1-calls/SAMPLE_${i}"
    export sample_name=$(cat "sample_name.txt")
    export container_name=${i}_${sample_name}
    export output_genotyped_intervals="/path/postprocess_case/${sample_name}_genotyped_intervals.vcf"
    export output_genotyped_segments="/path/postprocess_case/${sample_name}_genotyped_segments.vcf"
    export output_denoised_copy_ratios="/path/postprocess_case/${sample_name}_denoised_ratios.tsv"

    echo "#Sample index: $i, sample name: $sample_name" >> $run_file
    echo "docker run --name=$container_name -dt broadinstitute/gatk" >> $run_file

	echo "docker exec ${container_name} gatk PostprocessGermlineCNVCalls -R /path/reference.fna --contig-ploidy-calls /path/contig_ploidy_casemode/ploidy-case-calls --calls-shard-path /path/case_gCNV_calls/case_1-calls --model-shard-path /path/cohort_gCNV_calls/cohort_1-model ... --calls-shard-path /path/case_gCNV_calls/case_549-calls --model-shard-path /path/cohort_gCNV_calls/cohort_549-model --allosomal-contig chrX --allosomal-contig chrY --sample-index $i --output-genotyped-intervals $output_genotyped_intervals --output-genotyped-segments $output_genotyped_segments --output-denoised-copy-ratios $output_denoised_copy_ratios" >> $run_file
	echo "docker stop ${container_name}" >> $run_file
	echo "docker rm ${container_name}" >> $run_file
    echo "" >> $run_file
done

bash $run_file 2> /path/postprocessCNV_casemode_samples.log