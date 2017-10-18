#! /bin/bash 

# Automation script for running whole analysis. 
# bash WGH_Analysis <read_1.fq> <read_2.fq> <output_folder_name>

# Fetch paths to external software (defines WGH_path, bowtie2_ref, picard_path and fragScaff_path)
. paths.sh

# Argument parsing
r1=$1
r2=$2
output=$3

# Stripping extension from file name for adding others
file=$r1
name_ext=$(basename "$file")
name="${name_ext%.*}"
file_name="$path/${name_ext%.*}"

# Stripping extension from file name for adding others
file2=$r2
name_ext2=$(basename "$file2")
name2="${name_ext2%.*}"
file_name2="$path/${name_ext2%.*}"

#
# Intials
#

mkdir output

#
# Trimming & barcode extraction
#

# Trim away E handle on R1 5'.
cutadapt -g CAGTTGATCATCAGCAGGTAATCTGG \
    -o $file_name".h1_removed.fastq" \
    -p $file_name2".h1_removed.fastq" $r1 $r2 \
    --discard-untrimmed -e 0.2 &&

# Get DBS using UMI-Tools -> _BDHVBDVHBDVHBDVH in header.
umi_tools extract --stdin=$file_name".h1_removed.fastq" \
    --stdout=$file_name".cut_n_extract.fastq" \
    --bc-pattern=NNNNNNNNNNNNNNNNNNNN --bc-pattern2= \
    --read2-in=$file_name2".h1_removed.fastq" \
    --read2-out=$file_name2".cut_n_extract.fastq" \
    -L $file_name".cut_n_extract_log.txt" &&

#Cut TES from 5' of R1. TES=AGATGTGTATAAGAGACAG. Discard untrimmed.
cutadapt -g AGATGTGTATAAGAGACAG -o $file_name".TES1_removed.fastq" \
    -p $file_name2".TES1_removed.fastq" \
    $file_name".cut_n_extract.fastq" \
    $file_name2".cut_n_extract.fastq" --discard-untrimmed -e 0.2 &&

#Cut TES' from 3' for R1 and R2. TES'=CTGTCTCTTATACACATCT
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
	-o $file_name".TES2_removed.fastq" \
	-p $file_name2".TES2_removed.fastq" \
	$file_name".TES1_removed.fastq" \
	$file_name2".TES1_removed.fastq" -e 0.2 &&

#
# Mapping
#

bowtie2 --maxins 2000 -x /Users/pontushojer/data_analysis/Reference/Bowtie2Index/genome \
    -1 $r1 -2 $r2 -S $output/"mappedInserts.sam" | samtools view -bh - > $path/"mappedInserts.bam" &&

rm $path/"mappedInserts.sam" &&

samtools sort $path/"mappedInserts.bam" \
    -o $path/"mappedInserts.sort.bam"

#
# Clustering
#

python $WGH_path/python_scripts/cdhit_prep.py \
    -r 2 -f 0
for file in $output/indexing_barcodes/*.fa; do;
	cd-hit-454 -T 0 -c 0.9 -gap -100 -g 1 -n 3 -M 0 -i $file -o $file'.clustered';
	done
cat $output/indexing_barcodes/*.clstr > NN.clstr

#
# Tagging and duplcate removal
#

python $WGH_path/python_scripts/tag_bam.py 
picardtools rmdup

#
# Scaffolding
#

perl $fragScaff_path 
