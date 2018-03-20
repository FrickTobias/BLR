## Main steps
Here every step will be explained in detail, for exact commands used in-house see example folder. Whenever arrows (< >) 
encloses something, replace it with what is written inside (without the arrows). 

### Read trimming and barcode identification
This is comprised of four steps, steps 1-3 only trim read 1 and step 4 trims both sequences. First handle 1 is trimmed 
from the 5' end of read 1 followed by moving the barcode (20 bp) from the read sequence to the header. Subsequently 
h2 (TES) is trimmed from the same end and lastly, short inserts are trimmed for TES' in their 3' end. The handle 
sequences are as follows:

```
h1 = CAGTTGATCATCAGCAGGTAATCTGG
h2 = TES = CTGTCTCTTATACACATCT

```

```
# Trim handle 1
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
 
# Trim handle 2 
cutadapt -g AGATGTGTATAAGAGACAG \
    <read_1.h1.bc.fq> \    
    <read_2.h1.bc.fq> \
    -o <read_1.h1.bc.h2.fq> \
    -p <read_2.h1.bc.h2.fq> \
    -j <processors> \
    --discard_untrimmed \
    -g 0.2 \
    -m 65
    
# Trim handle 2:s reverse complement from 3' end (in both reads)  
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
    <read_1.h1.bc.h2.fq> \    
    <read_2.h1.bc.h2.fq> \
    -o <read_1.trimmed.fq> \
    -p <read_2.trimmed.fq> \
    -j <processors> \
    -e 0.2 \
    -m 25
```

### Mapping & Filtering

Here reads are mapped to the human genome. The mapping results are then filtered to remove unmapped reads and
non-primary alignements.

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
   
# Sorting file (necessary for downstream analysis)
samtools sort \
    <mapped.bam> > <mapped.sort.bam>
    
# Filtering for unmapped reads (0x04) and non-primary alignments (0xs100)
samtools view \
    <mapped.sort.bam> \
    -bh \
    -F 0x04 \
    -F 0x100 > <mapped.sort.filt.bam> 
    
```

### Clustering

The clustering is done in three steps; preparation, clustering and tagging. The preparation steps takes trimmed reads
and writes .fasta files containing only unique barcode sequences which then are clustered by CD-HIT-454. Lastly the
clustered barcodes are used are used to tag the .bam file with reads where they are given a cluster ID (a number) 
alongside their barcode sequence in the header as well as a RG tag with the barcode ID.

In order to cluster a vast amount of unique barcodes sequences (>100M) the barcode sequences are divided into several 
files according on its first nucleotides. These are used as an index and these files are clustered separately and this
can be reglated in the preparation step (cdhit_prep.py, -r option).

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

### Duplicate removal

Now the file is ready for removing read duplicates and cluster dupilcates. First read duplicates are removed (not if 
they have different barcodes) followed by marking remaining duplicates (not taking barcode into account), yielding a 
file with only cluster duplicates marked. The clusters marked as duplicates are then merged if two phased read pairs 
within phasing distance and they both share the same barcodes.

```
# Removing duplicates, keeping reads if they have different barcodes
java -jar picard.jar MarkDuplicates 
    I=<mapped.sort.filt.tag.bam>
    O=<mapped.sort.filt.tag.rmdup.bam>
    M=<rmdup.log>
    ASSUME_SORT_ORDER=true \
    REMOVE_DUPLICATES=true \
    BARCODE_TAG=RG
 
# Marking cluster duplicates
java -jar picard.jar MarkDuplicates
    I=<mapped.sort.filt.tag.rmdup.bam>
    O=<mapped.sort.filt.tag.rmdup.mkdup.bam>
    M=<mkdup.log>
   
# Merging cluster duplicates
python cluster_rmdup.py \
    <mapped.sort.filt.tag.rmdup.mkdup.bam> \
    <mapped.sort.filt.tag.rmdup.x2.bam>
```

By now you will have an analysis-ready file.