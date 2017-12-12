#!/bin/bash

#
# Initials
#

processors=1
mailing=false
remove=false

#
# Argument parsing
#

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
	    echo "  -r  removes files generated during analysis instead of just compressing them"
	    echo ''
	    exit 0
	    ;;
    esac
done

#
# Positonal redundancy for option useage
#

ARG1=${@:$OPTIND:1}
ARG2=${@:$OPTIND+1:1}

#
# Mailing option
#

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

#
# Error handling
#

if [ -z "$ARG1" ] || [ -z "$ARG2" ]
then
    echo ""
    echo "ARGUMENT ERROR"
    echo "Did not find all three positional arguments, see -h for more information."
    echo "(got bam:"$ARG1" and output:"$ARG2" instead)"
    echo ""
    exit 0
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

path=$ARG3

tagged_bam=ARG1
output=ARG2

java -jar $picard_path MarkDuplicates I=$bamfile O=$bamfile'.rmdup.bam' M=rmdup_log.txt ASSUME_SORT_ORDER=coordinate REMOVE_DUPLICATES=true BARCODE_TAG=RG

java -jar $picard_path MarkDuplicates I=$bamfile'.rmdup.bam' O=$bamfile'.rmdup.mkdup.bam' M=mkdup_cluster_log.txt ASSUME_SORT_ORDER=coordinate

python3 $wgh_path/python_scripts/cluster_rmdup.py $bamfile'.rmdup.mkdup.bam' $bamfile'.rmdup.x2.bam'