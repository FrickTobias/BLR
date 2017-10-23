#! /bin/bash

# Creates paths.txt file which can be run in another shell using '. paths.sh' to fetch the variables defined in the paths.sh script.

# bash setpath.txt <wgh_path> <picard_path> <bowtie2_reference> <fragScaff_path>

printf 'WGH_path='$1 > paths.txt
printf '\npicard_path='$2 >> paths.txt
printf '\nbowtie2_reference='$3 >> paths.txt
printf '\nfragScaff_path='$4 >> paths.txt
