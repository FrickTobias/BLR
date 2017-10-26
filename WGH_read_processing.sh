#!/bin/bash
#
# Script to process WFA reads without using the DBS_Analysis pipline
#

# RUN std
#    $ bash wfa_processing.sh <r1.fastq> <r2.fastq> <out.dir>
#
# RUN with log 
#    $ bash wfa_processing.sh <r1.fastq> <r2.fastq> <out.dir> > log.txt  2>$1

#
# Fetching paths to external programs (from paths.txt)
#

# PATH to WGH_Analysis folder
wgh_path=$(dirname "$0")

# Loading PATH:s to software
#   - reference:            $bowtie2_reference
#   - Picard tools:         $picard_path
#   - fragScaff:            $fragScafff_path
. $wgh_path'/paths.txt'

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

mkdir -p $path

#echo 'Starting 1st trim '$(date) | mail -s 'wgh' tobias.frick@scilifelab.se &&
printf '#1 START PROCESSING \n'

# Trim away E handle on R1 5'. Also removes reads shorter than 85 bp.
cutadapt -g ^CAGTTGATCATCAGCAGGTAATCTGG \
    -o $file_name".h1.fastq" \
    -p $file_name2".h1.fastq" $1 $2 \
    --discard-untrimmed -e 0.2 -m 65 # Tosses reads shorter than len(e+bc+handle+TES)

#echo 'Starting umi extraction '$(date) | mail -s 'wgh' tobias.frick@scilifelab.se &&
printf '\n\n#2 TRIMMED E \n'

# Get DBS using UMI-Tools -> _BDHVBDVHBDVHBDVH in header.
umi_tools extract --stdin=$file_name".h1.fastq" \
    --stdout=$file_name".h1.bc.fastq" \
    --bc-pattern=NNNNNNNNNNNNNNNNNNNN --bc-pattern2= \
    --read2-in=$file_name2".h1.fastq" \
    --read2-out=$file_name2".h1.bc.fastq" \
    -L $file_name".h1.bc.txt"

rm $file_name".h1.fastq"
rm $file_name2".h1.fastq"

#echo 'Starting 2nd trim '$(date) | mail -s 'wgh' tobias.frick@scilifelab.se &&
printf '\n\n#3 GOT DBS USING UMI-TOOLs \n'

#Cut TES from 5' of R1. TES=AGATGTGTATAAGAGACAG. Discard untrimmed.
cutadapt -g AGATGTGTATAAGAGACAG -o $file_name".h1.bc.h2.fastq" \
    -p $file_name2".h1.bc.h2.fastq" \
    $file_name".h1.bc.fastq" \
    $file_name2".h1.bc.fastq" --discard-untrimmed -e 0.2

rm $file_name".h1.bc.fastq"
rm $file_name2".h1.bc.fastq"

#echo 'Starting 3rd trim (final) '$(date) | mail -s 'wgh' tobias.frick@scilifelab.se &&
printf '\n\n#4 TRIMMED TES1 \n'

#Cut TES' from 3' for R1 and R2. TES'=CTGTCTCTTATACACATCT
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
	-o $file_name".trimmed.fastq" \
	-p $file_name2".trimmed.fastq" \
	$file_name".h1.bc.h2.fastq" \
	$file_name2".h1.bc.h2.fastq" -e 0.2


rm $file_name".h1.bc.h2.fastq"
rm $file_name2".h1.bc.h2.fastq"

#echo 'Starting mapping '$(date) | mail -s 'wgh' tobias.frick@scilifelab.se &&
printf '\n\n#5 TRIMMED TES2 \n'

bowtie2 --maxins 2000 -x $bowtie2_reference \
    -1 $file_name".trimmed.fastq" -2 $file_name2".trimmed.fastq" -S $path/"mappedInserts.sam"

samtools view -bh $path/"mappedInserts.sam" > $path/"mappedInserts.bam"

pigz $file_name".trimmed.fastq"
pigz $file_name2".trimmed.fastq"

#echo 'Starting bam file sorting '$(date) | mail -s 'wgh' tobias.frick@scilifelab.se &&
printf '\n\n#6 MAPPED READS \n'

rm $path/"mappedInserts.sam"

printf '\n\n#7 REMOVED SAM-FILE \n'

samtools sort $path/"mappedInserts.bam" \
    -o $path/"mappedInserts.sort.bam"

rm $path/"mappedInserts.bam"

samtools index $path/"mappedInserts.sort.bam"

#echo 'Finished '$(date) | mail -s 'wgh' tobias.frick@scilifelab.se &&
printf '\n\n#8 SORTED AND INDEXED BAM-FILE \n'

printf 'RUN COMPLETE (>'-')  (>'-')>  ^('-')^'
