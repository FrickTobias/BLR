#!/bin/bash

 #                                   #
############# Overview ################
 #                                   #
 # 1. argument parsing & initials    #
 #                                   #
 # 2. Start of script                #
 #   - Cut first handle (=e)         #
 #   - Extract barcode to header     #
 #   - Cut second handle (=TES)      #
 #   - Trim 3' end for TES'          #
 #                                   #
#######################################
 #                                   #

#
# Argument parsing & preparations
#

# Initials
processors=1
mailing=false
remove=false

# Parse arguments
while getopts "m:hp:r" OPTION
do
    case ${OPTION} in
    
        p)
            processors=${OPTARG}
            ;;
        m)
            email=${OPTARG}
            mailing=true
            ;;
        r)
            remove=true
            ;;
        h)
            printf '\nThis script runs the trimming parts of the WGH pipeline. Input are two WGH read files and output is written to a directory containing four sets of compressed fastq files. The final files are the ".trimmed.fq" files.

Useage: bash WGH_read_processing.sh <options> <r1.fq> <r2.fq> <output_dir>

Positional arguments (required)
  <r1.fq>         Read one in .fastq format. Also handles gzip files (.fastq.gz)
  <r2.fq>         Read two in .fastq format. Also handles gzip files (.fastq.gz)
  <output_dir>    Output directory for analysis results

Optional arguments
  -h  help (this output)
  -m  mails the supplied email when analysis is finished
  -p  processors for threading. DEFAULT: 1
  -r  removes files generated during analysis instead of just compressing them. DEFAULT: False

NB: options must be given before arguments.\n'
	        exit 0
	        ;;
    esac
done

# Positonal redundancy for option useage
ARG1=${@:$OPTIND:1}
ARG2=${@:$OPTIND+1:1}
ARG3=${@:$OPTIND+2:1}

# Mailing option
if $mailing
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

# Error handling
if [ -z "$ARG1" ] || [ -z "$ARG2" ] || [ -z "$ARG3" ]
then
    echo ""
    echo "ARGUMENT ERROR"
    echo "Did not find all three positional arguments, see -h for more information."
    echo "(got r1:"$ARG1", r2:"$ARG2" and output:"$ARG3" instead)"
    echo ""
    exit 0
fi

# Fetching paths to external programs (from paths.txt)

# PATH to WGH_Analysis folder
wgh_path=$(dirname "$0")

# Loading PATH:s to software
#   - reference:            $bowtie2_reference
#   - Picard tools:         $picard_path
. $wgh_path'/paths.txt'

# output folder
path=$ARG3

# File one prep
file=$ARG1
name_ext=$(basename "$file")
name="${name_ext%.*}"
file_name="$path/${name_ext%.*}"

# File two prep
file2=$ARG2
name_ext2=$(basename "$file2")
name2="${name_ext2%.*}"
file_name2="$path/${name_ext2%.*}"

# Mailing
if $mailing
    then
    echo 'Starting 1st trim '$(date) | mail -s 'wgh' $email
fi
printf '\nRunning with '$processors' threads\n'
printf '\n#1 START PROCESSING \n'

#
# Start of script
#

mkdir -p $path

# Trim away E handle on R1 5'. Also removes reads shorter than 85 bp.
cutadapt -g ^CAGTTGATCATCAGCAGGTAATCTGG \
    -j $processors \
    -o $file_name".h1.fastq" \
    -p $file_name2".h1.fastq" \
    $ARG1 \
    $ARG2 \
    --discard-untrimmed -e 0.2 -m 65 > $path/trimming.txt # Tosses reads shorter than len(e+bc+handle+TES)

printf '#2 TRIMMED E \n'
pigz $file_name".h1.fastq"
pigz $file_name2".h1.fastq"

# Get DBS using UMI-Tools -> _BDHVBDVHBDVHBDVH in header.
umi_tools extract --stdin=$file_name".h1.fastq.gz" \
    --stdout=$file_name".h1.bc.fastq" \
    --bc-pattern=NNNNNNNNNNNNNNNNNNNN --bc-pattern2= \
    --read2-in=$file_name2".h1.fastq.gz" \
    --read2-out=$file_name2".h1.bc.fastq" \
    -L $file_name".h1.bc.txt"

# Remove
if $remove
then
    rm $file_name".h1.fastq.gz"
    rm $file_name2".h1.fastq.gz"
fi

# Compress
pigz $file_name".h1.bc.fastq"
pigz $file_name2".h1.bc.fastq"
printf '#3 GOT DBS USING UMI-TOOLs \n'

#Cut TES from 5' of R1. TES=AGATGTGTATAAGAGACAG. Discard untrimmed.
cutadapt -g AGATGTGTATAAGAGACAG -o $file_name".h1.bc.h2.fastq" \
    -j $processors \
    -p $file_name2".h1.bc.h2.fastq" \
    $file_name".h1.bc.fastq.gz" \
    $file_name2".h1.bc.fastq.gz" \
    --discard-untrimmed -e 0.2  >> $path/trimming.txt

# Remove
if $remove
then
    rm $file_name".h1.bc.fastq.gz"
    rm $file_name2".h1.bc.fastq.gz"
fi

# Compress
pigz $file_name".h1.bc.h2.fastq"
pigz $file_name2".h1.bc.h2.fastq"
printf '#4 TRIMMED TES1 \n'

#Cut TES' from 3' for R1 and R2. TES'=CTGTCTCTTATACACATCT
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
    -j $processors \
	-o $file_name".trimmed.fastq" \
	-p $file_name2".trimmed.fastq" \
	-m 25 \
	$file_name".h1.bc.h2.fastq.gz" \
	$file_name2".h1.bc.h2.fastq.gz" \
	-e 0.2  >> $path/trimming.log

# Remove
if $remove
then
    rm $file_name".h1.bc.h2.fastq.gz"
    rm $file_name2".h1.bc.h2.fastq.gz"
fi

# Compress/remove
pigz $file_name".trimmed.fastq"
pigz $file_name2".trimmed.fastq"

if $mailing
    then
    echo 'Trimming finished '$(date) | mail -s $path $email
fi

# Ugly solution to calculate % construOK
var1=$( cat $path/trimming.txt | grep 'Read 1 with adapter' | cut -d '(' -f 2 | cut -d '%' -f 1 | tr '\n' ' ' | cut -d ' ' -f 1 )
var2=$( cat $path/trimming.txt | grep 'Read 1 with adapter' | cut -d '(' -f 2 | cut -d '%' -f 1 | tr '\n' ' ' | cut -d ' ' -f 2 )

printf 'RUN COMPLETE\n'
awk '{print "\nIntact reads: "$1*$2*0.0001" %\n"}' <<< "$var1 $var2"
