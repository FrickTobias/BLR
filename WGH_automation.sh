#! /bin/bash 

# Automation script for running whole analysis. 
# bash WGH_Analysis <read_1.fq> <read_2.fq> <output_folder_name>

# Fetch paths to external software (defines WGH_path, bowtie2_ref, picard_path and fragScaff_path)
. paths.sh

# Argument parsing
r1=$1
r2=$2
output=$3

#
# Trimming
#

bash $WGH_path/bash_scripts/cutadapt1.sh 
bash $WGH_path/bash_scripts/UMItools_extract.sh
bash $WGH_path/bash_scripts/cutadapt2.sh

#
# Clustering
#

python $WGH_path/python_scripts/cdhit_prep.py -r 2 -f 0 
for file in $output/indexing_barcodes/*.fa; do;
	cd-hit-454 -T 0 -c 0.9 -gap -100 -g 1 -n 3 -M 0 -i $file -o $file'.clustered';
	done
cat $output/indexing_barcodes/*.clstr > NN.clstr

#
# Mapping
#

bowtie2 

#
# Tagging and duplcate removal
#

python $WGH_path/python_scripts/tag_bam.py 
picardtools rmdup

#
# Scaffolding
#

perl $fragScaff_path 
