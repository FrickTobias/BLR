#!/bin/bash

#
# Argument parsing & preparations
#

# Initials
processors=1
mailing=false
remove=false
heapspace=90

# Parse arguments
while getopts "m:hp:r:H" OPTION
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
        H)
            heapspace=${OPTARG}
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
  -H  maximum heapspace (~RAM) for duplicate removal step. DEFAULT: 90Gb
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
    echo 'Emailing '$email ' when finished'
    echo 'Starting run '$(date) | mail -s 'wgh' $email
fi

#
# Start of script
#

printf '\n Starting analysis'
printf 'MAKE SURE YOU HAVE AT LEAST 6X THE ORGINAL FILE SIZE AVAILABLE DISK SPACE'
printf 'MAKE SURE YOUR SYSTEM HAS '$heapspace ' GB RAM AVAILABLE'


. wgh_path/WGH_read_processing.sh -r -p $processors $ARG1 $ARG2 $path

. wgh_path/WGH_read_mapper.sh -r -p $processors $path/$file_name".trimmed.fastq.gz" $path/$file_name2".trimmed.fastq.gz" $path/map.filt.sort.bam

. wgh_path/WGH_cluster_bc.sh -r -p $processors $path/$file_name".trimmed.fastq.gz" $path/map.filt.sort.bam $path/unique_bc

. wgh_path/WGH_rmdup.sh -r -H $heapspace'G' -p $processors $path/map.filt.sort.bam.tag.bam $path/map.filt.sort.bam.tag.rmdup.x2.bam

if $mailing
then
    echo 'Analysis finished '$date | mail -s 'wgh' $email
fi