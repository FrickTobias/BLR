# Creates paths.txt file which can be run in another shell using '. paths.sh' to fetch the variables defined in the paths.sh script.

# bash setpath.txt <picard_path> <bowtie2_reference> <fragScaff_path>

echo 'YOU HAVE TO BE IN THE WGH_ANALYSIS FOLDER FOR THIS TO WORK!'
echo ''
echo 'This will be patched shortly'
echo 'As of now, control that paths.txt has the correct stucture manually'

printf '\npicard_path='$1 > paths.txt
printf '\nbowtie2_reference='$2 >> paths.txt
printf '\nfragScaff_path='$3 >> paths.txt
