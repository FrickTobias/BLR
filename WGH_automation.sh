#!/bin/bash

# WGH_automation
#
# author: Tobias Frick
# mail: tobias.frick@scilifelab.se
# github: https://github.com/FrickTobias/WGH_Analysis

  #                                         #
# # # # # # # # # # # # # # # # # # # # # # # #
  #                                         #
  # 0. Argument parsing & initials          #
  # 1. Read trimming & demultiplexing       #
  # 2. Mapping & filtering                  #
  # 3. Barcode clustering & tagging         #
  # 4. Clstr rmdup & fastq generation       #
  #   - optional, requires A LOT of RAM     #
  #                                         #
# # # # # # # # # # # # # # # # # # # # # # # #
  #                                         #


# 0. ###################################################################################

 #                                   #
############# Overview ################
 #                                   #
 #   - Argument parsing              #
 #   - Option handling               #
 #   - Variable names & paths        #
 #                                   #
#######################################
 #                                   #

# Initials
processors=1
mailing=false
remove=false
duplicate_rmdup=false
keep_logiles_true=false

# Argparsing
while getopts "m:hp:r:d:k" OPTION
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
        d)
            duplicate_rmdup=true
            ;;
        k)
            keep_logfiles=true
            ;;
        h)
            printf '\nThis script runs the trimming parts of the WGH pipeline. Input are two WGH read files and output is written to a directory containing four sets of compressed fastq files. The final files are the ".trimmed.fq" files.

Useage: bash WGH_read_processing.sh <options> <r1.fq> <r2.fq> <output_dir>
Example: bash WGH_read_proccessing.sh -k -r -p 24 -m john.doe@myemail.com human_N1298_read1.fastq human_N1298_read2.fastq human_N1298

Positional arguments (required)
  <r1.fq>         Read one in .fastq format. Also handles gzip files (.fastq.gz)
  <r2.fq>         Read two in .fastq format. Also handles gzip files (.fastq.gz)
  <output_dir>    Output directory for analysis results

Optional arguments
  -h  help (this output)
  -m  mails the supplied email when analysis is finished
  -p  processors for threading. DEFAULT: 1
  -r  removes files generated during analysis instead of just compressing them. DEFAULT: False
  -d  duplicates (read & cluster) removal step will also be run. Not suited for RAM < 100GB. DEFAULT: False
  -k  keep all logfiles generated during analysis instead of keeping only specifics. DEFAULT: False

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

# Processor option
printf '\nRunning with '$processors' threads\n'

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
mkdir -p $path

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

# Logfiles
trim_logfile = $path'/1_trim.log'
map_logfile = $path'/2_map.log'
cluster_logfile = $path'/3_cluster.log'
if duplicate_rmdup
then
    rmdup_logfile = $path'/4_rmdup.log'
fi

# Mailing
if $mailing
then
    echo 'ANALYSIS STARTING '$(date) | mail -s $path $email
fi

# 1. ###################################################################################

 #                                   #
############# Overview ################
 #                                   #
 #   - Cut first handle (=e)         #
 #   - Extract barcode to header     #
 #   - Cut second handle (=TES)      #
 #   - Trim 3' end for TES'          #
 #                                   #
#######################################
 #                                   #

# Mailing
if $mailing
then
    echo '1_trim starting '$(date) | mail -s $path $email
fi
printf "`date`"'\nRemoving 1st adaptor\n'

# Trim away E handle on R1 5'. Also removes reads shorter than 85 bp.
cutadapt -g ^CAGTTGATCATCAGCAGGTAATCTGG \
    -j $processors \
    -o $file_name".h1.fastq" \
    -p $file_name2".h1.fastq" \
    $ARG1 \
    $ARG2 \
    --discard-untrimmed -e 0.2 -m 65 > $trim_logfile # Tosses reads shorter than len(e+bc+handle+TES)

printf "`date`"'\n1st adaptor removal done\n'
printf "`date`"'\nExtracting barcode\n'

## Get DBS using UMI-Tools -> _BDHVBDVHBDVHBDVH in header.
python3 $wgh_path'/python scripts/bc_extract.py' \
    $file_name".h1.fastq"
    $file_name2".h1.fastq"
    $file_name".h1.bc.fastq"
    $file_name2"h1.bc.fastq"
if $remove
then
    rm $file_name".h1.fastq.gz"
    rm $file_name2".h1.fastq.gz"
fi
pigz $file_name".h1.fastq"
pigz $file_name2".h1.fastq"
pigz $file_name".h1.bc.fastq"
pigz $file_name2".h1.bc.fastq"

printf "`date`"'\nBarcode extraction done\n'
printf "`date`"'\nRemoving 2nd adaptor\n'

#Cut TES from 5' of R1. TES=AGATGTGTATAAGAGACAG. Discard untrimmed.
cutadapt -g AGATGTGTATAAGAGACAG -o $file_name".h1.bc.h2.fastq" \
    -j $processors \
    -p $file_name2".h1.bc.h2.fastq" \
    $file_name".h1.bc.fastq.gz" \
    $file_name2".h1.bc.fastq.gz" \
    --discard-untrimmed -e 0.2  >> $trim_logfile
if $remove
then
    rm $file_name".h1.bc.fastq.gz"
    rm $file_name2".h1.bc.fastq.gz"
fi
pigz $file_name".h1.bc.h2.fastq"
pigz $file_name2".h1.bc.h2.fastq"

printf "`date`"'\n2nd adaptor removed\n'
printf "`date`""\nTrimming 5'"

#Cut TES' from 3' for R1 and R2. TES'=CTGTCTCTTATACACATCT
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
    -j $processors \
	-o $file_name".trimmed.fastq" \
	-p $file_name2".trimmed.fastq" \
	-m 25 \
	$file_name".h1.bc.h2.fastq.gz" \
	$file_name2".h1.bc.h2.fastq.gz" \
	-e 0.2  >> $path/trimming.log
if $remove
then
    rm $file_name".h1.bc.h2.fastq.gz"
    rm $file_name2".h1.bc.h2.fastq.gz"
fi
pigz $file_name".trimmed.fastq"
pigz $file_name2".trimmed.fastq"

# Ugly solution to calculate % construOK
var1=$( cat $trim_logfile | grep 'Read 1 with adapter' | cut -d '(' -f 2 | cut -d '%' -f 1 | tr '\n' ' ' | cut -d ' ' -f 1 )
var2=$( cat $trim_logfile | grep 'Read 1 with adapter' | cut -d '(' -f 2 | cut -d '%' -f 1 | tr '\n' ' ' | cut -d ' ' -f 2 )
awk '{print "\nIntact reads: "$1*$2*0.0001" %\n"}' <<< "$var1 $var2"

if $mailing
    then
    echo '1_trim finished '$(date) | mail -s $path $email
fi
printf "`date`""\n5' trimmed"

# 2. ###################################################################################

 #                                   #
############# Overview ################
 #                                   #
 #   - Map & convert to bam          #
 #   - Sorting                       #
 #   - Filtering (unmap + prim map)  #
 #                                   #
#######################################
 #                                   #

if $mailing
then
    echo '2_map starting '$(date) | mail -s $path $email
fi

printf "`date`"'\tMapping starting\n'
pritnf '\n\n Map stats: .sort.bam\n' >> $map_logfile

# Mapping & bam conversion
(bowtie2 \
    -1 $file_name".trimmed.fastq.gz" \
    -2 $file_name2".trimmed.fastq.gz" \
    -x $bowtie2_reference \
    --maxins 2000 \
    -p $processors | \
    samtools view \
        - \
        -@ $processors \
        -bh > $file_name".bam") 2>$map_logfile

printf "`date`"'\tMapping done\n'
printf "`date`"'\tSorting\n'

# Sorting
samtools sort \
    $file_name".bam" \
    -@ processors > $file_name".sort.bam"

if $remove
then
    rm $file_name".bam"
fi

printf "`date`"'\tSorting done\n'
printf "`date`"'\tMap stats\n'
pritnf '\n\n flagstat: .bam\n' >> $map_logfile
samtools flagstat \
    $file_name".sort.bam" >> $map_logfile


printf "`date`"'\tMap stats done\n'
printf "`date`"'\tFiltering\n'

# Filtering
samtools view \
    $file_name".sort.bam" \
    -@ $processors \
    -bh \
    -F 0x04 \
    -F 0x100 > $file_name".sort.filt.bam"

printf "`date`"'\tFiltering done\n'
printf "`date`"'\tMap stats\n'

pritnf '\n\n flagstat: sort.filt.bam\n' >> $map_logfile
samtools flagstat \
    $file_name".sort.filt.bam" >> $map_logfile

if $mailing
then
    echo '2_map finished '$(date) | mail -s $path $email
fi

printf "`date`"'\tMap stats done\n'

# 3. ###################################################################################

 #                                   #
############# Overview ################
 #                                   #
 #   - Bc seq extraction             #
 #   - Clustering                    #
 #   - Merge clust files & tag bam   #
 #                                   #
#######################################
 #                                   #

if $mailing
then
    echo '3_clustering starting'$(date) | mail -s $path $email
fi
printf "`date`"'\tExtracting barcodes\n'

# Barcode extraction
pigz -d $file_name".trimmed.fastq.gz"
python3 $wgh_path'/python scripts/cdhit_prep.py' \
    $file_name".trimmed.fastq" \
    $path"unique_bc" \
    -r 3\
    -f 0
pigz $file_name".trimmed.fastq"

printf "`date`"'\tBarcode extraction done\n'
printf "`date`"'\tClustering barcodes\n'

# Barcode clustering
for file in $path"/unique_bc"/*.fa
do
    printf '\n' >> $map_logfile
    wc -l $file >> $map_logfile
    cd-hit-454 \
        -i $file \
        -o $file'.clustered' \
        -T $processors \
        -c 0.9 \
        -gap 100 \
        -g 1 \
        -n 3 \
        -M 0
done
cat $path"/unique_bc/"*".clstr" > $path"/NNN.clstr"

if $remove
then
    rm -rf $path"/unique_bc"
fi

printf "`date`"'\tBarcodes clustering done\n'
printf "`date`"'\tTaggin bam\n'

# Tagging bamfile
python3 $wgh_path'/python scripts/tag_bam.py' \
    $file_name".sort.filt.bam" \
    $path/NNN.clstr \
    $file_name".sort.filt.tag.bam"

if $mailing
then
    echo '3_clustering finished '$(date) | mail -s $path $email
fi
printf "`date`"'\tBam taggin done\n'

# 4. ###################################################################################

 #                                   #
############# Overview ################
 #                                   #
 #   - rmdup (with TAG=RG)           #
 #   - mkdup (without TAG)           #
 #   - Cluster rmdup                 #
 #   - Fastq generation              #
 #                                   #
#######################################
 #                                   #

if duplicate_rmdup:
then

    printf "`date`"'\tRemoving duplicates\n'

    if $mailing
    then
        echo '4_rmdup starting '$(date) | mail -s $path $email
    fi


    java -jar $picard_path MarkDuplicates \
        I=$file_name".sort.filt.tag.bam" \
        O=$file_name".sort.filt.tag.rmdup.bam" \
        M=$rmdup_logfile \
        ASSUME_SORT_ORDER=coordinate \
        REMOVE_DUPLICATES=true BARCODE_TAG=RG

    printf "`date`"'\tDuplicate removal done\n'
    printf "`date`"'\tMarking barcode duplicates\n'

    java -jar $picard_path MarkDuplicates \
        I=$file_name".sort.filt.tag.rmdup.bam" \
        O=$file_name".sort.filt.tag.rmdup.mkdup.bam" \
        M=$path/picard_mkdup.txt \
        ASSUME_SORT_ORDER=coordinate
    cat $path/picard_mkdup.txt >> $rmdup_logfile
    rm $path/picard_mkdup.txt

    printf "`date`"'\tBarcode duplicates marking done\n'
    printf "`date`"'\tMerging clusters\n'

    python3 $wgh_path'/python scripts/cluster_rmdup.py' \
        $file_name".sort.filt.tag.rmdup.mkdup.bam" \
        $file_name".sort.filt.tag.rmdup.x2.bam"

    printf "`date`"'\tCluster merging done\n'
    printf "`date`"'\tGenerating fastqs\n'

    java -jar $picard_path SamToFastq \
        I=$file_name".sort.filt.rmdup.x2.bam" \
        FASTQ=$file_name".final.fastq" \
        SECOND_END_FASTQ=$file_name2".final.fastq"

    pigz $file_name".final.fastq"
    pgiz $file_name2".final.fastq"

    if $mailing
    then
        echo '4_rmdup finished '$(date) | mail -s $path $email
    fi
    printf "`date`"'\tFastq generation done\n'

fi

printf "`date`"'\tFINISHED\n'

if $mailing
then
    echo 'ANALYSIS FINISHED '$(date) | mail -s $path $email
fi