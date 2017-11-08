#!/bin/bash

#
# Initials
#

processors=1
mailing=False
remove=False

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
            mailing=True
            ;;
        r)
            remove=True
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

path=$ARG3

file=$ARG1
name_ext=$(basename "$file")
name="${name_ext%.*}"
file_name="$path/${name_ext%.*}"

file2=$ARG2
name_ext2=$(basename "$file2")
name2="${name_ext2%.*}"
file_name2="$path/${name_ext2%.*}"

mkdir -p $path

# Mailing
if $mailing
    then
    echo 'Starting 1st trim '$(date) | mail -s 'wgh' $email
fi
printf '#1 START PROCESSING \n'

# Trim away E handle on R1 5'. Also removes reads shorter than 85 bp.
cutadapt -g ^CAGTTGATCATCAGCAGGTAATCTGG \
    -o $file_name".h1.fastq" \
    -p $file_name2".h1.fastq" $ARG1 $ARG2 \
    --discard-untrimmed -e 0.2 -m 65 # Tosses reads shorter than len(e+bc+handle+TES)

# Mailing
if $mailing
    then
    echo 'Starting umi extraction '$(date) | mail -s 'wgh' $email
fi
printf '\n\n#2 TRIMMED E \n'

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

# Mailing
if $mailing
    then
    echo 'Starting 2nd trim '$(date) | mail -s 'wgh' $email
fi

printf '\n\n#3 GOT DBS USING UMI-TOOLs \n'

#Cut TES from 5' of R1. TES=AGATGTGTATAAGAGACAG. Discard untrimmed.
cutadapt -g AGATGTGTATAAGAGACAG -o $file_name".h1.bc.h2.fastq" \
    -p $file_name2".h1.bc.h2.fastq" \
    $file_name".h1.bc.fastq.gz" \
    $file_name2".h1.bc.fastq.gz" --discard-untrimmed -e 0.2

# Remove
if $remove
then
    rm $file_name".h1.bc.fastq.gz"
    rm $file_name2".h1.bc.fastq.gz"
fi

# Compress
pigz $file_name".h1.bc.h2.fastq"
pigz $file_name2".h1.bc.h2.fastq"


# Mailing
if $mailing
    then
    echo 'Starting 3rd trim (final) '$(date) | mail -s 'wgh' $email
fi

printf '\n\n#4 TRIMMED TES1 \n'

#Cut TES' from 3' for R1 and R2. TES'=CTGTCTCTTATACACATCT
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
	-o $file_name".trimmed.fastq" \
	-p $file_name2".trimmed.fastq" \
	-m 25 \
	$file_name".h1.bc.h2.fastq.gz" \
	$file_name2".h1.bc.h2.fastq.gz" -e 0.2

# Remove
if $remove
then
    rm $file_name".h1.bc.h2.fastq.gz"
    rm $file_name2".h1.bc.h2.fastq.gz"
fi

# Compress/remove
pigz $file_name".trimmed.fastq"
pigz $file_name2".trimmed.fastq"

#echo 'Starting mapping '$(date) | mail -s 'wgh' tobias.frick@scilifelab.se
#printf '\n\n#5 TRIMMED TES2 \n'
#
#bowtie2 --maxins 2000 -x $bowtie2_reference \
#    -1 $file_name".trimmed.fastq" -2 $file_name2".trimmed.fastq" | samtools view -bS - > $path/"mappedInserts.bam"
#
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

if $mailing
    then
    echo 'Finished '$(date) | mail -s 'wgh' $email
fi

#printf '\n\n#8 SORTED AND INDEXED BAM-FILE \n'
#
#printf 'RUN COMPLETE (>'-')  (>'-')>  ^('-')^'
