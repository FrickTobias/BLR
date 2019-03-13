## Main steps
Here every step will be explained in detail with descriptions of commands, options and functionality for all pipeline
steps. Whenever arrows (< >) encloses something, replace it with what is written inside (without the arrows). 

### Read trimming and barcode identification
This is comprised of four steps, steps 1-3 trim read 1 and step 4 trims both read 1 and read 2. First handle 1 is trimmed 
from the 5' end of read 1 followed by moving the barcode (20 bp) from the read sequence to the header. Subsequently 
handle 2 is trimmed from the same end and lastly, sequences are trimmed for the reverse complement of the handles 
flanking the insert gDNA from their 3'-ends. This is to ensure read structures on the other end is not present, in case
of short insert gDNA sequences (less than 215). Handle sequences are as follows:

```
handle 1 = CAGTTGATCATCAGCAGGTAATCTGG
handle 2 = CTGTCTCTTATACACATCT

```

And is trimmed using these four commands:

```
# Trim handle 1 and discard untrimmed read pairs.
cutadapt -g ^CAGTTGATCATCAGCAGGTAATCTGG \
    <read_1.fq> \
    <read_2.fq> \
    -o <read_1.h1.fq> \
    -p <read_2.h1.fq> \
    -j <processors> \
    --discard_untrimmed \
    -e 0.2 \
    -m 65 

# Remove and extract bc to header  
python3 bc_extract.py \
    <read_1.h1.fq> \
    <read_2.h1.fq> \
    <read_1.h1.bc.fq> \
    <read_2.h1.bc.fq>
 
# Trim handle 2 and discard untrimmed read pairs (keeps untrimmed)
cutadapt -g AGATGTGTATAAGAGACAG \
    <read_1.h1.bc.fq> \    
    <read_2.h1.bc.fq> \
    -o <read_1.h1.bc.h2.fq> \
    -p <read_2.h1.bc.h2.fq> \
    -j <processors> \
    --discard_untrimmed \
    -e 0.2 \
    
# Trim handle 2:s reverse complement from 3' end in both reads  
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
    <read_1.h1.bc.h2.fq> \    
    <read_2.h1.bc.h2.fq> \
    -o <read_1.trimmed.fq> \
    -p <read_2.trimmed.fq> \
    -j <processors> \
    -e 0.2 \
    -m 25
```

Options used:

```
   OPTION           FUNCTION
   -e 0.2           Demands at least 20% similarity to sequence handles for positive match
   -m 65            Demands a minimum of 65 remaining bp after trimming
   -g ^H1           Looks for H1 at the start (^) of 5' in read 1 (g)
```

### Mapping & Filtering

Here reads are mapped to the human genome and converted from sam to bam format.

```
# Mapping reads and writing as bam
bowtie2 \
    -1 <read_1.trimmed.fq> \
    -2 <read_2.trimmed.fq> \
    -x <reference/Bowtie2/genome> \
    --maxins 2000 \
    -p <processors> |
    samtools view \
        - \
        -@ <processors>
        -bh > <mapped.bam>
   
# Sorting file (necessary for downstream analysis, will also reduce memory needed)
samtools sort \
    <mapped.bam> > <mapped.sort.bam>
        
```

Options used:

```
   OPTION           FUNCTION
   --maxins 2000    Maximum distance read_1<->read_2 can cover.
   -b               Writes as .bam (compressed version of the standard format .sam)
   -h               Include header in output, default is to not write it
```

### Clustering

The clustering is done in three steps; preparation, clustering and tagging. The preparation steps takes trimmed reads
and writes fasta files containing unique barcode sequences which are clustered separately by CD-HIT-454. Lastly 
the clustered barcodes are used are used to tag mapped reads with a barcode ID (corresponding to clustering results) 
alongside their barcode sequence in the header as well as a BC tag with the barcode ID.

In order to cluster a vast amount of unique barcodes sequences (>100M), barcode sequences are divided into several 
files according on its first nucleotides. This can be regulated by using the cdhit_prep.py -r (--reduce) option or
in the -i (--indexing_nucleotides) option in the automation script.

```
# Preparation, using 3 indexing nucleotides
python cdhit_prep.py \
    <read_1.fq> \
    <output_barcode_directory>
    -r <index>
    -f 0
    
# Clustering all indexing files separately
for index_file in <output_barode_directory>/*.fa
do
    cd-hit-454 \
        -i $index_file \
        -o $index_file'.clustered' \
        -T <processors> \
        -c 0.9 \
        -gap 100 \
        -g 1 \
        -n 3 \
        -M 0
done

# Merging results into one file
cat <outpute_barcode_directory>/*.clstr > barcodes_<index>.clstr
 
# Tagging .bam file with clustered barcodes
python3 tag_bam.py \
    mapped.sort.filt.bam \
    barcodes_<index>.clstr \
    mapped.sort.filt.tag.bam

```
Options used:

```
   OPTION           FUNCTION
   -r 3             Uses first three bases to divide barcodes into separate files.
   -f 0             Do not use a read count threshold to include barcode sequences.
   -gap 100         High gap opening score
   -g 1             The more accurate but slower method
   -n 3             k-mer lenth for similarity definition
   -M 0             Unlimited memory
```

### Duplicate removal

Now the file is ready for removing read duplicates and cluster duplicates. First read duplicates are removed (not if 
they have different barcodes) followed by marking remaining duplicates (not taking barcode into account), yielding a 
file with only read duplicates with different barcodes marked. These reads are then used to find pairs of phased reads
within the same clusters which have another pair of phased reads at the same position with another barcode (pairs of 
phased reads are defined as two proximal read pairs). If this is identified this is interpreted as the barcodes either 
originating from the same emulsion or not being missed in the clustering step. Following this step the file can have 
its barcodes stripped for reads in droplets with more than -M molecules. 

```
# Removing duplicates, keeping reads if they have different barcodes
java -jar picard.jar MarkDuplicates 
    I=<mapped.sort.filt.tag.bam>
    O=<mapped.sort.filt.tag.rmdup.bam>
    M=<rmdup.log>
    ASSUME_SORT_ORDER=true \
    REMOVE_DUPLICATES=true \
    BARCODE_TAG=BC
 
# Marking cluster duplicates
java -jar picard.jar MarkDuplicates
    I=<mapped.sort.filt.tag.rmdup.bam>
    O=<mapped.sort.filt.tag.rmdup.mkdup.bam>
    M=<mkdup.log>
   
# Merging cluster duplicates
python cluster_rmdup.py \
    <mapped.sort.filt.tag.rmdup.mkdup.bam> \
    <mapped.sort.filt.tag.rmdup.x2.bam>

# Filter big clusters
python filter_clusters.py \
    <mapped.sort.filt.tag.rmdup.x2.bam> \
    <stats_prefix> \
    -f <mapped.sort.filt.tag.rmdup.x2.filt.bam> \
    -M 260
```
Options used:

```
   OPTION           FUNCTION
   ASO=coordinate   Assume file is sorted based on mapping positions
   BARCODE_TAG=BC   Bamfile tag used for storing barcoding information
   -f <file>        Write an output bam file with big droplets filtered (barcodes with many molecules)
   -M 260           Maximum amount of molecules allowed per barcode
   
```

By now you will have an analysis-ready file.