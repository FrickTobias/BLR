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
# Initials
#

processors=1
mailing=False

#
# Argument parsing
#

while getopts "m:hp:" OPTION; do
	case ${OPTION} in

	    p)
	        processors={OPTARG}
		    ;;
        m)
            email={OPTARG}
            mailing=True
            ;;
		h)
		    echo ''
			echo 'This script runs the trimming parts of the WGH pipeline. Input are two WGH read files and output is written to a directory containing four sets of compressed fastq files. The final files are the ".trimmed.fq" files.'
			echo ""
			echo 'Useage: bash WGH_read_processing.sh <r1.fq> <r2.fq> <output_dir>'
			echo ''
			echo "Positional arguments (required)"
			echo "  <r1.fq>         Read one in .fastq format. Also handles gzip files (.fastq.gz)"
			echo "  <r2.fq>         Read two in .fastq format. Also handles gzip files (.fastq.gz)"
			echo "  <output_dir>    Output directory for analysis results"
			echo ""
			echo "Optional arguments"
			echo "  -h  help (this output)"
			echo "  -m  mails the supplied email when analysis is finished"
			echo "  -p  processors for threading, not implemented yet"
			echo ''
			exit 0
			;;
	esac
done

if [ $mailing == True ]
        then
        if [[ $email == *"@"* ]]
                then
                echo 'Mailing '$email' when finished.'
        else
                echo 'FORMAT ERROR: -m '
                echo ''
                echo 'Please supply email on format john.doe@domain.org'
                echo '(got "'$email'" instead)'
                exit 0
        fi
fi

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

mkdir -p $path

if [ $mailing == True ]
    then
    echo 'Starting 1st trim '$(date) | mail -s 'wgh' $email
fi
printf '#1 START PROCESSING \n'

# Trim away E handle on R1 5'. Also removes reads shorter than 85 bp.
cutadapt -g ^CAGTTGATCATCAGCAGGTAATCTGG \
    -o $file_name".h1.fastq" \
    -p $file_name2".h1.fastq" $1 $2 \
    --discard-untrimmed -e 0.2 -m 65 # Tosses reads shorter than len(e+bc+handle+TES)

if [ $mailing == True ]
    then
    echo 'Starting umi extraction '$(date) | mail -s 'wgh' $email
fi
printf '\n\n#2 TRIMMED E \n'

# Get DBS using UMI-Tools -> _BDHVBDVHBDVHBDVH in header.
umi_tools extract --stdin=$file_name".h1.fastq" \
    --stdout=$file_name".h1.bc.fastq" \
    --bc-pattern=NNNNNNNNNNNNNNNNNNNN --bc-pattern2= \
    --read2-in=$file_name2".h1.fastq" \
    --read2-out=$file_name2".h1.bc.fastq" \
    -L $file_name".h1.bc.txt"

pigz $file_name".h1.fastq"
pigz $file_name2".h1.fastq"

if [ $mailing == True ]
    then
    echo 'Starting 2nd trim '$(date) | mail -s 'wgh' $email
fi

printf '\n\n#3 GOT DBS USING UMI-TOOLs \n'

#Cut TES from 5' of R1. TES=AGATGTGTATAAGAGACAG. Discard untrimmed.
cutadapt -g AGATGTGTATAAGAGACAG -o $file_name".h1.bc.h2.fastq" \
    -p $file_name2".h1.bc.h2.fastq" \
    $file_name".h1.bc.fastq" \
    $file_name2".h1.bc.fastq" --discard-untrimmed -e 0.2

pigz $file_name".h1.bc.fastq"
pigz $file_name2".h1.bc.fastq"

if [ $mailing == True ]
    then
    echo 'Starting 3rd trim (final) '$(date) | mail -s 'wgh' $email
fi
printf '\n\n#4 TRIMMED TES1 \n'

#Cut TES' from 3' for R1 and R2. TES'=CTGTCTCTTATACACATCT
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
	-o $file_name".trimmed.fastq" \
	-p $file_name2".trimmed.fastq" \
	-m 25 \
	$file_name".h1.bc.h2.fastq" \
	$file_name2".h1.bc.h2.fastq" -e 0.2

pigz $file_name".h1.bc.h2.fastq"
pigz $file_name2".h1.bc.h2.fastq"

#echo 'Starting mapping '$(date) | mail -s 'wgh' tobias.frick@scilifelab.se
#printf '\n\n#5 TRIMMED TES2 \n'
#
#bowtie2 --maxins 2000 -x $bowtie2_reference \
#    -1 $file_name".trimmed.fastq" -2 $file_name2".trimmed.fastq" -S $path/"mappedInserts.sam"
#
#samtools view -bh $path/"mappedInserts.sam" > $path/"mappedInserts.bam"
#
#pigz $file_name".trimmed.fastq"
#pigz $file_name2".trimmed.fastq"
#
#echo 'Starting bam file sorting '$(date) | mail -s 'wgh' tobias.frick@scilifelab.se
#printf '\n\n#6 MAPPED READS \n'
#
#rm $path/"mappedInserts.sam"
#
#printf '\n\n#7 REMOVED SAM-FILE \n'
#
#samtools sort $path/"mappedInserts.bam" \
#    -o $path/"mappedInserts.sort.bam"
#
#rm $path/"mappedInserts.bam"
#
#samtools index $path/"mappedInserts.sort.bam"
#

if [ $mailing == True ]
    then
    echo 'Finished '$(date) | mail -s 'wgh' $email
fi

#printf '\n\n#8 SORTED AND INDEXED BAM-FILE \n'
#
#printf 'RUN COMPLETE (>'-')  (>'-')>  ^('-')^'
