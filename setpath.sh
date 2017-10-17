#! /bin/bash

# Creates paths.sh file which can be run in another shell using '. paths.sh' to fetch the variables defined in the paths.sh script.

# bash setpath.sh <picard_path> <bowtie2_reference> <fragScaff_path>

printf '#! /bin/bash/\n' > paths.sh
printf '\npicard_path='$1'' >> paths.sh
printf '\nbowtie2_reference='$2 >> paths.sh
printf '\nfragScaff_path='$3 >> paths.sh
