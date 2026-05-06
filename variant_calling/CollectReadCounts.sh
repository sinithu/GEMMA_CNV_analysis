#!/usr/bin/env bash

# https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-rare-germline-copy-number-variants
# Counting read starts for all 1000 bp intervals
# Performs some filtering, such as removal of duplicate reads and reads with mapping quality under 10

# This script first outputs a run file with CollectReadCounts commands for all the samples, and then executes it

export run_file="/path/CollectReadCounts_samples.sh"
echo "#!/usr/bin/env bash" > $run_file
sample_files=$(ls /path/bamfiles/*.bam)

for file in $sample_files
do
    name=$(basename $file .bam)
    echo "#$name" >> $run_file
    echo "docker run --name=$name -dt broadinstitute/gatk" >> $run_file
    echo "docker exec $name gatk CollectReadCounts -L /path/preprocessed_intervals.interval_list -imr OVERLAPPING_ONLY -I $file -O /path/counts/${name}.hdf5" >> $run_file
    echo "docker stop $name" >> $run_file
	echo "docker rm $name" >> $run_file
    echo "" >> $run_file
done

bash $run_file 2> /path/CollectReadCounts_samples.log