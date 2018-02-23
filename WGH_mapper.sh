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
	    echo "  <output_bam>    Bamfile with mapped reads."
	    echo ""
	    echo "Optional arguments"
	    echo "  -h  help (this output)"
	    echo "  -m  mails the supplied email when analysis is finished"
	    echo "  -p  processors for threading. DEFAULT: 1"
	    echo "  -r  removes files generated during analysis"
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
ARG3=${@:$OPTIND+2:1}

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

if [ -z "$ARG1" ] || [ -z "$ARG2" ] || [ -z "$ARG3" ]
then
    echo ""
    echo "ARGUMENT ERROR"
    echo "Did not find all three positional arguments, see -h for more information."
    echo "(got r1:"$ARG1", r2:"$ARG2" and output:"$ARG3" instead)"
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

r1=$ARG1
r2=$ARG2
output_bam=$ARG3

if $mailing
then
    echo 'starting mapping '$date | mail -s 'wgh' $email
fi

printf "`date`"'\tMapping starting\n'

(bowtie2 --maxins 2000 -p $processors -x $bowtie2_reference \
    -1 $r1 \
    -2 $r2 | \
    samtools view -@ $processors -bh - > $output_bam.tmp.bam) 2>mapper.log

printf "`date`"'\tMapping done\n'
printf "`date`"'\tFiltering starting\n'

printf '\n'mapped'\n' >> mapper.log
samtools flagstat $output_bam.tmp.bam >> mapper.log

samtools view -@ $processors -bh -F 0x04 -F 0x100 $output_bam.tmp.bam > $output_bam.tmp.filt.bam

if $remove
then
    rm $output_bam.tmp.bam
fi

printf "`date`"'\tFiltering done\n'
printf "`date`"'\tSorting starting\n'

printf '\n'mapped.filt'\n' >> mapper.log
samtools flagstat $output_bam.tmp.filt.bam >> mapper.log

samtools sort -@ processors $output_bam.tmp.filt.bam > $output_bam

if $remove
then
    rm $output_bam.tmp.filt.bam
fi

printf "`date`"'\tSorting done\n'
printf "`date`"'\tFINISHED\n'


if $mailing
then
    echo 'mapping finished '$date | mail -s 'wgh' $email
fi
