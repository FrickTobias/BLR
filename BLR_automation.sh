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
keep_logiles=false
heap_space=90
index_nucleotides=3
threshold=0
start_step=1
end_step=4

# Argparsing
while getopts "hrkh:m:p:i:H:s:e:t:" OPTION
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
        k)
            keep_logfiles=true
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
  -m  mails the supplied email when analysis is finished                                DEFAULT: None
  -p  processors for threading                                                          DEFAULT: 1
  -r  removes files generated during analysis instead of just compressing them          DEFAULT: false
  -h  help (this output)                                                                DEFAULT: N/A

Advanced options: globals
  -s  start at this step number (see Pipeline outline)                                  DEFAULT: 1
  -e  end after this step number (see Pipeline outline)                                 DEFAULT: 4
  -k  keep all logfiles generated during analysis instead of keeping only specifics     DEFAULT: false

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


# Mailing option
if $mailing
then
    if [[ $email == *"@"* ]]
    then
        printf '\nMail:\t\t'$email
    else
        echo ''
        echo 'OPTION ERROR: -m '
        echo ''
        echo 'Please supply email on format john.doe@domain.org'
        echo '(got "'$email'" instead)'
        exit 0
    fi
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

# Mailing
if $mailing
then
    echo 'ANALYSIS STARTING '$(date) | mail -s $path $email
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

    # Mailing
    if $mailing
    then
        echo '1_trim starting '$(date) | mail -s $path $email
    fi
    printf '\n1. Demultiplexing\n'
    printf "`date`"'\t1st adaptor removal\n'

    # Trim away E handle on R1 5'. Also removes reads shorter than 85 bp.
    cutadapt -g ^CAGTTGATCATCAGCAGGTAATCTGG \
        -j $processors \
        -o $file_name".h1.fastq.gz" \
        -p $file_name2".h1.fastq.gz" \
        $ARG1 \
        $ARG2 \
        --discard-untrimmed -e 0.2 -m 65 > $trim_logfile # Tosses reads shorter than len(e+bc+handle+TES)

    printf "`date`"'\t1st adaptor removal done\n'
    printf "`date`"'\tBarcode extraction\n'

    ## Get DBS using UMI-Tools -> _BDHVBDVHBDVHBDVH in header.
    (python3 $wgh_path'/python scripts/bc_extract.py' \
        $file_name".h1.fastq.gz" \
        $file_name2".h1.fastq.gz" \
        $file_name".h1.bc.fastq" \
        $file_name2".h1.bc.fastq") 2>$path"/bc_extract.stderr"
    if ! $keep_logiles
    then
        rm $path"/bc_extract.stderr"
    fi
    if $remove
    then
        rm $file_name".h1.fastq.gz"
        rm $file_name2".h1.fastq.gz"
    fi
    pigz $file_name".h1.bc.fastq"
    pigz $file_name2".h1.bc.fastq"

    printf "`date`"'\tBarcode extraction done\n'
    printf "`date`"'\t2nd adaptor removal\n'

    #Cut TES from 5' of R1. TES=AGATGTGTATAAGAGACAG. Discard untrimmed.
    cutadapt -g AGATGTGTATAAGAGACAG \
        -o $file_name".h1.bc.h2.fastq.gz" \
        -j $processors \
        -p $file_name2".h1.bc.h2.fastq.gz" \
        $file_name".h1.bc.fastq.gz" \
        $file_name2".h1.bc.fastq.gz" \
        --discard-untrimmed -e 0.2  >> $trim_logfile
    if $remove
    then
        rm $file_name".h1.bc.fastq.gz"
        rm $file_name2".h1.bc.fastq.gz"
    fi

    printf "`date`"'\t2nd adaptor removal done\n'
    printf "`date`""\t3' trimming\n"

    #Cut TES' from 3' for R1 and R2. TES'=CTGTCTCTTATACACATCT
    cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
        -j $processors \
        -o $file_name".trimmed.fastq.gz" \
        -p $file_name2".trimmed.fastq.gz" \
        -m 25 \
        $file_name".h1.bc.h2.fastq.gz" \
        $file_name2".h1.bc.h2.fastq.gz" \
        -e 0.2  >> $trim_logfile


    if $remove
    then
        rm $file_name".h1.bc.h2.fastq.gz"
        rm $file_name2".h1.bc.h2.fastq.gz"
    fi

    if $mailing
    then
        echo '1_trim finished '$(date) | mail -s $path $email
    fi
    printf "`date`""\t3' trimming done\n"

    # Ugly solution to calculate % construOK
    var1=$( cat $trim_logfile | grep 'Read 1 with adapter' | cut -d '(' -f 2 | cut -d '%' -f 1 | tr '\n' ' ' | cut -d ' ' -f 1 )
    var2=$( cat $trim_logfile | grep 'Read 1 with adapter' | cut -d '(' -f 2 | cut -d '%' -f 1 | tr '\n' ' ' | cut -d ' ' -f 2 )
    printf "`date`""\t"; awk '{print "Intact reads: "$1*$2*0.0001" %"}' <<< "$var1 $var2"

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

    if $mailing
    then
        echo '2_clustering starting'$(date) | mail -s $path $email
    fi
    printf '\n2. Clustering\n'
    printf "`date`"'\tBarcode fasta generation\n'

    # Barcode extraction
    pigz -d $file_name".trimmed.fastq.gz"
    (python3 $wgh_path'/python scripts/cdhit_prep.py' \
        $file_name".trimmed.fastq" \
        $path"/unique_bc" \
        -i $index_nucleotides\
        -f 0 >$path"/cdhit_prep.stdout") 2>$path"/cdhit_prep.stderr"
    if ! $keep_logiles
    then
        rm $path"/cdhit_prep.stdout"
        rm $path"/cdhit_prep.stderr"
    fi
    pigz $file_name".trimmed.fastq"

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

    if ! $keep_logiles
    then
        rm $path"/cdhit.log"
    fi

    if $remove
    then
        rm -rf $path"/unique_bc"
    fi

    if $mailing
    then
        echo '2_clustering finished '$(date) | mail -s $path $email
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

    if $mailing
    then
        echo '3_map starting '$(date) | mail -s $path $email
    fi
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
    (python3 $wgh_path'/python scripts/tag_bam.py' \
        $file_name".sort.bam" \
        $path"/"$N_string".clstr" \
        $file_name".sort.tag.bam" ) 2>$path"/tag_bam.stderr"

    if ! $keep_logiles
    then
        rm $path"/tag_bam.stderr"
    fi

    if $mailing
    then
        echo '3_map finished '$(date) | mail -s $path $email
    fi

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

    if $mailing
    then
        echo '4_rmdup starting '$(date) | mail -s $path $email
    fi
    printf '\n4. Duplicate removal\n'
    printf "`date`"'\tDuplicate removal\n'

    # Read duplicate removal
    (java '-Xmx'$heap_space'g' -jar $picard_path MarkDuplicates \
        I=$file_name".sort.tag.bam" \
        O=$file_name".sort.tag.rmdup.bam" \
        M=$path"/picard.log" \
        ASSUME_SORT_ORDER=coordinate \
        REMOVE_DUPLICATES=true \
        BARCODE_TAG=BC) 2>$rmdup_logfile

    printf "`date`"'\tDuplicate removal done\n'
    printf "`date`"'\tBarcode duplicate marking\n'

    # Cluster duplicate marking
    (java '-Xmx'$heap_space'g' -jar $picard_path MarkDuplicates \
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
    (java -jar $picard_path SamToFastq \
        I=$file_name".sort.tag.rmdup.x2.filt.bam" \
        FASTQ=$file_name".final.fastq" \
        SECOND_END_FASTQ=$file_name2".final.fastq") 2>>$path/cpicard.log

    if ! $keep_logiles
    then
        rm $file_name".sort.tag.rmdup.x2.bam.log"
        rm $path/picard.log
    fi

    pigz $file_name".final.fastq"
    pigz $file_name2".final.fastq"

    if $mailing
    then
        echo '4_rmdup finished '$(date) | mail -s $path $email
    fi
    printf "`date`"'\tFastq generation done\n'

fi

printf '\n'"`date`"'\tANALYSIS FINISHED\n'

if $mailing
then
    echo 'ANALYSIS FINISHED '$(date) | mail -s $path $email
fi