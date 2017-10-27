#!/usr/bin/env bash

# PATH to WGH_Analysis folder
wgh_path=$(dirname "$0")

# Loading PATH:s to software
#   - reference:            $bowtie2_reference
#   - Picard tools:         $picard_path
#   - fragScaff:            $fragScafff_path
. $wgh_path'/paths.txt'

r1 = $1
bam = $2
output = $3

echo 'Starting bc extraction '$(date) | mail -s 'wgh' tobias.frick@scilifelab.se &&
python3 $wgh_path/python_scripts/cdhit_prep.py $r1 $output/unique_bc -r 3 -f 0

echo 'Starting clustering '$(date) | mail -s 'wgh' tobias.frick@scilifelab.se &&
for file in $output/unique_bc/*;
    do;
    cdhit-454 -i $file -o $file'.clustered' -T 0 -c 0.9 -gap -100 -g 1 -n 3 -M 0;
    done;

echo 'Starting tagging '$(date) | mail -s 'wgh' tobias.frick@scilifelab.se &&
cat $output/unique_bc/*.clstr > $output/NNN.clstr

python3 $wgh_path/python_scripts/tag_bam.py $output/NNN.clstr $bam $bam'.tagged.bam'

echo 'Finished '$(date) | mail -s 'wgh' tobias.frick@scilifelab.se &&