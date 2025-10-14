#!/usr/bin/env bash

# https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-rare-germline-copy-number-variants
# For COHORT MODE
# Calling CNVs across scatters in all samples
# This script first outputs a run file with GermlineCNVCaller commands for all the scatters, and then executes it

#First, creating a variable with all paths of scatters from /path/scatter
#For 549 scatters
scatters=$(for i in {0001..0549}
do
	echo /path/scatter/temp_${i}_of_549/scattered.interval_list
done)

#Creating a GermlineCNVCaller file for all scatters and samples
#All read count files from CollectReadCounts need to be added with "-I" flags
export run_file="/path/GermlineCNVCaller_cohortmode_scatters.sh"
echo "#!/usr/bin/env bash" > $run_file
export count=0

for scatter in $scatters
do
    count=$(($count+1))
    export container_name="gCNV_cohort_${count}"

    echo "#Scatter number $count" >> $run_file
    echo "docker run --name=$container_name -dt broadinstitute/gatk" >> $run_file
    echo "docker exec $container_name gatk GermlineCNVCaller --run-mode COHORT -L $scatter --contig-ploidy-calls /path/contig_ploidy/ploidy-calls --annotated-intervals /path/annotated_intervals.tsv -imr OVERLAPPING_ONLY --output /path/cohort_gCNV_calls --output-prefix cohort_${count} --verbosity DEBUG -I /path/counts/sample1.hdf5 -I /path/counts/sample2.hdf5 -I ... -I /path/counts/sampleN.hdf5" >> $run_file
    echo "docker stop $container_name" >> $run_file
	echo "docker rm $container_name" >> $run_file
    echo "" >> $run_file
done

bash $run_file 2> /path/GermlineCNVCaller_cohortmode_scatters.log