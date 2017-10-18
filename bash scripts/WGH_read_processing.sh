#!/bin/bash
#
# Script to process WFA reads without using the DBS_Analysis pipline
#

# RUN std
#    $ bash wfa_processing.sh <r1.fastq> <r2.fastq> <out.dir>
#
# RUN with log 
#    $ bash wfa_processing.sh <r1.fastq> <r2.fastq> <out.dir> > log.txt  2>$1

path=$3 

file=$1
name_ext=$(basename "$file")
name="${name_ext%.*}"
file_name="$path/${name_ext%.*}"

file2=$2
name_ext2=$(basename "$file2")
name2="${name_ext2%.*}"
file_name2="$path/${name_ext2%.*}"

#log="$path/log_file.txt"

mkdir -p $path &&

printf '#1 START PROCESSING \n' &&


# Trim away E handle on R1 5'.
cutadapt -g CAGTTGATCATCAGCAGGTAATCTGG \
    -o $file_name".e_removed.fastq" \
    -p $file_name2".e_removed.fastq" $1 $2 \
    --discard-untrimmed -e 0.2 &&


printf '\n\n#2 TRIMMED E \n' &&


# Get DBS using UMI-Tools -> _BDHVBDVHBDVHBDVH in header.
umi_tools extract --stdin=$file_name".e_removed.fastq" \
    --stdout=$file_name".cut_n_extract.fastq" \
    --bc-pattern=NNNNNNNNNNNNNNNNNNNN --bc-pattern2= \
    --read2-in=$file_name2".e_removed.fastq" \
    --read2-out=$file_name2".cut_n_extract.fastq" \
    -L $file_name".cut_n_extract_log.txt" &&


printf '\n\n#3 GOT DBS USING UMI-TOOLs \n' &&


#Cut TES from 5' of R1. TES=AGATGTGTATAAGAGACAG. Discard untrimmed.
cutadapt -g AGATGTGTATAAGAGACAG -o $file_name".TES1_removed.fastq" \
    -p $file_name2".TES1_removed.fastq" \
    $file_name".cut_n_extract.fastq" \
    $file_name2".cut_n_extract.fastq" --discard-untrimmed -e 0.2 &&


printf '\n\n#4 TRIMMED TES1 \n' &&


#Cut TES' from 3' for R1 and R2. TES'=CTGTCTCTTATACACATCT
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
	-o $file_name".TES2_removed.fastq" \
	-p $file_name2".TES2_removed.fastq" \
	$file_name".TES1_removed.fastq" \
	$file_name2".TES1_removed.fastq" -e 0.2 &&

printf '\n\n#5 TRIMMED TES2 \n' &&

bowtie2 --maxins 2000 -x /Users/pontushojer/data_analysis/Reference/Bowtie2Index/genome \
    -1 $1 -2 $2 -S $path/"mappedInserts.sam" | samtools view -bh - > $path/"mappedInserts.bam" &&

printf '\n\n#6 MAPPED READS \n' &&

rm $path/"mappedInserts.sam" &&

printf '\n\n#7 REMOVED SAM-FILE \n' &&

samtools sort $path/"mappedInserts.bam" \
    -o $path/"mappedInserts.sort.bam" &&

samtools index $path/"mappedInserts.sort.bam" &&

printf '\n\n#8 SORTED AND INDEXED BAM-FILE \n' &&

printf 'RUN COMPLETE (>'-')  (>'-')>  ^('-')^'
