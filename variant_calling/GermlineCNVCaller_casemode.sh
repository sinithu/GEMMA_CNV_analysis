#!/usr/bin/env bash

# https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-rare-germline-copy-number-variants
# For CASE MODE
# Calling CNVs across scatters in samples
# This script first outputs a run file with GermlineCNVCaller commands for all the scatters, and then executes it

#Creating a GermlineCNVCaller file for all scatters and samples
#All read count files from CollectReadCounts need to be added with "-I" flags
#For 549 scatters/models
export run_file="/path/GermlineCNVCaller_casemode_scatters.sh"
echo "#!/usr/bin/env bash" > $run_file

for i in {1..549}
do
    export container_name="gCNV_case_${i}"

    echo "#Scatter number $i" >> $run_file
    echo "docker run --name=$container_name -dt broadinstitute/gatk" >> $run_file
    echo "docker exec $container_name gatk GermlineCNVCaller --run-mode CASE --contig-ploidy-calls /path/contig_ploidy_case/ploidy-case-calls --model /path/cohort_gCNV_calls/cohort_${i}-model --output /path/case_gCNV_calls --output-prefix case_${i} --verbosity DEBUG -I /path/counts/sample1.hdf5 -I /path/counts/sample2.hdf5 -I ... -I /path/counts/sampleN.hdf5" >> $run_file
    echo "docker stop $container_name" >> $run_file
	echo "docker rm $container_name" >> $run_file
    echo "" >> $run_file
done

bash $run_file 2> /path/GermlineCNVCaller_casemode_scatters.log