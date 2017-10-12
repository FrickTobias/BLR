#!/bin/bash

# bash index_clustering.sh <output_directory> <read_1> <mapped_inserts>
# bash index_clustering.sh results r1_with_bc.fq mapped_inserts_with_bc.bam

python_folder=~/PycharmProjects/DBS_scratchbook/
bash_folder=~/bash_scripts/

output_directory=$1/
r1=$2
mapped_inserts=$3

python $python_folder/cdhit_prep.py -f 0 -r 2 $r1 $output_directory/unique_bc

for file in $output_directory/unique_bc/*.fa; do bash $bash_folder/cdhit.sh $file; done

cat $output_directory/unique_bc/*.clstr > output_directory/unique_bc/NN.clstr

python $python_folder/tag_bam.py $mapped_inserts $output_directory/unique_bc/NN.clstr $output_directory/$mapped_inserts'.tagged.bam'
