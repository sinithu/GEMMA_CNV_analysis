cd /path/filtered_VCFs

for file in *.gz
do
	name=$(basename $file _genotyped_segments_filtered.vcf.gz)
	docker run --name=$name -dt ensemblorg/ensembl-vep
	docker exec $name vep -i /path/filtered_VCFs/$file --format vcf --vcf -o /path/VCFs/${name}_genotyped_segments_filtered_VEP.vcf --fasta /path/reference.fna --offline --cache --everything --fork 4 --dir_cache /path/dir_cache
	docker stop $name
	docker rm $name
	bgzip /path/VCFs/${name}_genotyped_segments_filtered_VEP.vcf
done
