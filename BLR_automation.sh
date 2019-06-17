#!/bin/bash
set -euo pipefail

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
  # 4. Rmdup + filtering + fq-generation    #
  #                                         #
# # # # # # # # # # # # # # # # # # # # #na # # #
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
remove=false
duplicate_rmdup=false
heap_space=90
index_nucleotides=3
threshold=0
start_step=1
end_step=4

# Argparsing
while getopts "hrh:m:p:i:H:s:e:t:" OPTION
do
    case ${OPTION} in

        p)
            processors=${OPTARG}
            ;;
        r)
            remove=true
            ;;
        i)
            index_nucleotides=${OPTARG}
            ;;
        H)
            heap_space=${OPTARG}
            ;;
        t)
            threshold=${OPTARG}
            ;;
        s)
            start_step=${OPTARG}
            ;;
        e)
            end_step=${OPTARG}
            ;;
        h)
            printf 'BLR_automation.sh

Useage:     bash BLR_automation.sh <options> <r1.fq> <r2.fq> <output_dir>
Example:    bash BLR_automation.sh -r -p 24 -m john.doe@domain.org N1298_read1.fastq N1298_read2.fastq 180220_N1298
NB:         options must be given before arguments.

Pipeline outline:
  0.Argparsing & options
  1.Demultiplexing
  2.Clustering
  3.Mapping
  4.Duplicate removal

Positional arguments (REQUIRED)
  <r1.fq>       Read one in .fastq format. Also handles gzip files (.fastq.gz)
  <r2.fq>       Read two in .fastq format. Also handles gzip files (.fastq.gz)
  <output_dir>  Output directory for analysis results

Global optional arguments
  -p  processors for threading                                                          DEFAULT: 1
  -r  removes files generated during analysis instead of just compressing them          DEFAULT: false
  -h  help (this output)                                                                DEFAULT: N/A

Advanced options: globals
  -s  start at this step number (see Pipeline outline)                                  DEFAULT: 1
  -e  end after this step number (see Pipeline outline)                                 DEFAULT: 4

Advanced options: software settings
  -i  indexing nucletide number used for clustering (cdhit_prep.py)                     DEFAULT: 3
  -t  threshold for cluster duplicate calling (cluster_rmdup.py)                        DEFAULT: 0
  -H  heap space (~RAM) in GB for duplicate removal step                                DEFAULT: 90
  \n'
	        exit 0
	        ;;
    esac
done

# Positonal redundancy for option useage
ARG1=${@:$OPTIND:1}
ARG2=${@:$OPTIND+1:1}
ARG3=${@:$OPTIND+2:1}

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

printf '\n0. Argparsing & options'
printf '\nRead 1:\t'$ARG1'\nRead 2:\t'$ARG2'\nOutput:\t'$ARG3'\n'
printf '\nThreads:\t'$processors
printf '\nStarts at step:\t'$start_step
printf '\nEnd after step:\t'$end_step


# Fetching paths to external programs (from paths.txt)

# PATH to WGH_Analysis folder
wgh_path=$(dirname "$0")

# Loading PATH:s to software
#   - reference:            $bowtie2_reference
#   - Picard tools:         $picard_path
. $wgh_path'/paths.txt'

if [ -z ${picard_command+x} ]; then
    picard_command="java -Xmx${heap_space}g -jar $picard_path"
fi

# output folder
path=$PWD/$ARG3
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
trim_logfile=$path'/1_trim.log'
cluster_logfile=$path'/2_cluster.log'
map_logfile=$path'/3_map.log'
rmdup_logfile=$path'/4_rmdup.log'

# Remaining options
current_step=0
continue=true
if (( "$start_step" > "$end_step" ))
then
    printf "\n\nOPTION ERROR\nStart step cannot be larger than end step, see -s and -e option\n"
    exit 0
elif (( "$start_step" < 1 )) || (( "$start_step" > 4 ))
then
    printf "\n\nOPTION ERROR\nStart step must be within 1 and 4\n"
    exit 0
fi

# Make barcode file according BC.clstr OR BC.NNN.clstr, where N will correspond to how many index bases are used.
    if [[ $index_nucleotides == 0 ]]
    then
        N_string='BC'
    else
        N_string='BC.'
        for i in $(seq 1 $index_nucleotides)
        do
            N_string=$N_string'N'
        done
    fi

printf '\n\n'"`date`"'\tANALYSIS STARTING\n'

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

# Check if this step should be run
current_step=$((current_step+1))
if (( "$current_step" >= "$start_step" )) && [ "$continue" == true ]
then

    printf '\n1. Demultiplexing\n'
    printf "`date`"'\t1st adaptor removal\n'

    ln -s $PWD/$ARG1 $path/reads.1.fastq.gz
    ln -s $PWD/$ARG2 $path/reads.2.fastq.gz
    snakemake $path/trimmed-c.1.fastq.gz $path/trimmed-c.2.fastq.gz
    if $remove
    then
        rm $file_name".h1.bc.fastq.gz"
        rm $file_name2".h1.bc.fastq.gz"
    fi

    ln -s $path/trimmed-c.1.fastq.gz $file_name".trimmed.fastq.gz"
    ln -s $path/trimmed-c.2.fastq.gz $file_name2".trimmed.fastq.gz"

    printf "`date`""\t3' trimming done\n"

    # Ugly solution to calculate % construOK
    # TODO
#    var1=$( cat $trim_logfile | grep 'Read 1 with adapter' | cut -d '(' -f 2 | cut -d '%' -f 1 | tr '\n' ' ' | cut -d ' ' -f 1 )
#    var2=$( cat $trim_logfile | grep 'Read 1 with adapter' | cut -d '(' -f 2 | cut -d '%' -f 1 | tr '\n' ' ' | cut -d ' ' -f 2 )
#    printf "`date`""\t"; awk '{print "Intact reads: "$1*$2*0.0001" %"}' <<< "$var1 $var2"

fi

if (( "$current_step" == "$end_step" ))
then
    continue=false
fi


# 2. ###################################################################################

 #                                   #
############# Overview ################
 #                                   #
 #   - Bc seq extraction             #
 #   - Clustering                    #
 #   - Merge clust files & tag bam   #
 #                                   #
#######################################
 #                                   #

# Check if this step should be run
current_step=$((current_step+1))
if (( "$current_step" >= "$start_step" )) && [ "$continue" == true ]
then

    printf '\n2. Clustering\n'
    printf "`date`"'\tBarcode fasta generation\n'

    # Barcode extraction
    blr cdhitprep \
        $file_name".trimmed.fastq.gz" \
        $path"/unique_bc" \
        -i $index_nucleotides\
        -f 0 >$path"/cdhit_prep.stdout" 2>$path"/cdhit_prep.stderr"

    printf "`date`"'\tBarcode fasta generation done\n'
    printf "`date`"'\tBarcode clustering\n'

    # Non-indexing primer run fix
    if [[ $index_nucleotides == 0 ]]
    then
        mv $path"/unique_bc" $path"/unique_bc.fa"
        mkdir $path"/unique_bc"
        mv $path"/unique_bc.fa" $path"/unique_bc/unique_bc.fa"
    fi

    # Barcode clustering
    touch $path"/cdhit.log"
    for file in $path"/unique_bc"/*.fa
    do
        printf '\n' >> $cluster_logfile
        wc -l $file >> $cluster_logfile
        (cd-hit-454 \
            -i $file \
            -o $file'.clustered' \
            -T $processors \
            -c 0.9 \
            -gap 100 \
            -g 1 \
            -n 3 \
            -M 0) >> $path"/cdhit.log"
    done

    cat $path"/unique_bc/"*".clstr" > $path"/"$N_string".clstr"

    if $remove
    then
        rm -rf $path"/unique_bc"
    fi

    printf "`date`"'\tBarcode clustering done\n'

fi

if (( "$current_step" == "$end_step" ))
then
    continue=false
fi

# 3. ###################################################################################

 #                                   #
############# Overview ################
 #                                   #
 #   - Map & convert to bam          #
 #   - Sorting                       #
 #   - Filtering (unmap + prim map)  #
 #                                   #
#######################################
 #                                   #

# Check if this step should be run
current_step=$((current_step+1))
if (( "$current_step" >= "$start_step" )) && [ "$continue" == true ]
then

    printf '\n3. Mapping\n'
    printf "`date`"'\tMapping\n'
    printf '\n\n Map stats: .sort.bam\n' >> $map_logfile

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
    printf "`date`"'\tBam tagging\n'

    # Tagging bamfile
    (blr tagbam ' \
        $file_name".sort.bam" \
        $path"/"$N_string".clstr" \
        $file_name".sort.tag.bam" ) 2>$path"/tag_bam.stderr"

    printf "`date`"'\tBam tagging done\n'

fi

if (( "$current_step" == "$end_step" ))
then
    continue=false
fi

# 4. ###################################################################################

 #                                   #
############# Overview ################
 #                                   #
 #   - rmdup (with TAG=RG)           #
 #   - mkdup (without TAG)           #
 #   - Cluster rmdup                 #
 #   - Cluster filtering             #
 #   - Fastq generation              #
 #                                   #
#######################################
 #                                   #

# Check if this step should be run
current_step=$((current_step+1))
if (( "$current_step" >= "$start_step" )) && [ "$continue" == true ]
then

    printf '\n4. Duplicate removal\n'
    printf "`date`"'\tDuplicate removal\n'

    # Read duplicate removal
    ($picard_command MarkDuplicates \
        I=$file_name".sort.tag.bam" \
        O=$file_name".sort.tag.rmdup.bam" \
        M=$path"/picard.log" \
        ASSUME_SORT_ORDER=coordinate \
        REMOVE_DUPLICATES=true \
        BARCODE_TAG=BC) 2>$rmdup_logfile

    printf "`date`"'\tDuplicate removal done\n'
    printf "`date`"'\tBarcode duplicate marking\n'

    # Cluster duplicate marking
    ($picard_command MarkDuplicates \
        I=$file_name".sort.tag.rmdup.bam" \
        O=$file_name".sort.tag.rmdup.mkdup.bam" \
        M=$path"/mkdup.log" \
        ASSUME_SORT_ORDER=coordinate) 2>>$rmdup_logfile
    cat $path/"/mkdup.log" >> $path"/picard.log"
    rm $path"/mkdup.log"

    if $remove
    then
        rm $file_name".sort.tag.rmdup.bam"
    fi

    printf "`date`"'\tBarcode duplicate marking done\n'
    printf "`date`"'\tCluster merging\n'

    # Cluster duplicate merging
    (python3 $wgh_path'/python scripts/cluster_rmdup.py' \
        $file_name".sort.tag.rmdup.mkdup.bam" \
        $file_name".sort.tag.rmdup.x2.bam") 2>>$rmdup_logfile

    printf "`date`"'\tCluster merging done\n'
    printf "`date`"'\tIndexing\n'

    samtools index $file_name".sort.tag.rmdup.x2.bam"

    printf "`date`"'\tIndexing done\n'
    printf "`date`"'\tCluster filtering\n'

    mkdir -p $path"/cluster_stats"
    # Cluster filtering
    (python3 $wgh_path'/python scripts/filter_clusters.py' \
        -f $file_name".sort.tag.rmdup.x2.filt.bam" \
        -M 260 \
        $file_name".sort.tag.rmdup.x2.bam" \
        $path"/cluster_stats/x2.stats") 2>>$rmdup_logfile

    printf "`date`"'\tCluster filtering done\n'
    printf "`date`"'\tFastq generation\n'

    # Fastq generation
    ($picard_command SamToFastq \
        I=$file_name".sort.tag.rmdup.x2.filt.bam" \
        FASTQ=$file_name".final.fastq" \
        VALIDATION_STRINGENCY=SILENT \
        SECOND_END_FASTQ=$file_name2".final.fastq") 2>>$path/cpicard.log

    pigz $file_name".final.fastq"
    pigz $file_name2".final.fastq"

    printf "`date`"'\tFastq generation done\n'

fi

printf '\n'"`date`"'\tANALYSIS FINISHED\n'
