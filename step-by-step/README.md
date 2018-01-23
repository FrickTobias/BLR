### Main steps
Here, every step will be explained in detail. For more details about in-house option useage, see example folder. There
are also semi-automated scripts for performing the main parts of the analysis, see [semi_auto.md](https://github.com/FrickTobias/WGH_Analysis/blob/master/step-by-step/semi_automated.md).

#### Read trimming and barcode identification
This is comprised of four steps. First handle 1 is trimmed from the 5' end of read1. Then the next 20 bases are moved 
from the sequence up to the header using UMItools extract followed by another trimming of another hadle.

```
cutadapt \
    [-OPTIONS] \
    <read_1.fq> \
    <read_2.fq> \
    <read_1.h1.fq> \
    <read_2.h1.fq> 
    
UMItools \
    [-OPTIONS] \
    <read_1.h1.fq> \
    <read_2.h1.fq> \
    <read_1.h1.bc.fq> \
    <read_2.h1.bc.fq>
    
cutadapt \
    [-OPTIONS] \ 
    <read_1.h1.bc.fq> \    
    <read_2.h1.bc.fq> \
    <read_1.h1.bc.h2.fq> \
    <read_2.h1.bc.h2.fq>
    
cutadapt \
    [-OPTIONS] \ 
    <read_1.h1.bc.h2.fq> \    
    <read_2.h1.bc.h2.fq> \
    <read_1.trimmed.fq> \
    <read_2.trimmed.fq> 
```

#### Mapping

Mapping is done with bowtie2.

```
bowtie2 \
    [-OPTIONS] \
    <read_1.fq> \
    <read_2.fq> \
    <reference/Bowtie2/genome> \
    <output>
```

#### Clustering
The clustering is based on dividing barcode sequences in several files according to their first bases. First, 
barcode sequences are extracted from one of the read file headers and divided into a number of files (-r), 
depending on how many indexing bases are used. Optimally none would be used but using more significantly increases 
clustering speed. After this, run CD-HIT-454 on all generated fasta files in the generated directory and merge clstr
files to one single file.

Option recommendation for cdhit_prep.py:

   - 2 base indexing (-r 2)

```
python cdhit_prep.py \
    [-OPTIONS] \
    <read_1.fq> \
    <output_barcode_directory>
    
for index_file in <barcode_directory>/*.clstr
do
    cd-hit-454 7
        [-OPTIONS] \
        -i $index_file \
        -o $index_file'.clustered' 
done
    
cat <barcode>/.clstr > barcodes.clstr
```

#### Filtering

The mapped files are first filtered for mapping flags (unmapped and non-primary alignments) using samtools followed by 
being sorted. Afterwards duplicates are removed but without removing duplicates tagged with different barcode
ID:s (found in the @RG tag). 

Following this, a necessary preparation step is to mark duplicates but this time without 
but without taking barcode ID:s intpo account. Now the file will be ready for cluster duplicatate
removal which cluster_rmdup.py is used for.

```
samtools view
    [-OPTIONS]
    <inserts.tag.bam> 
    > <inserts.tag.filt.bam>
    
samtools sort \
    [-OPTIONS]
    <inserts.tag.filt.bam> \
    > <inserts.tag.filt.sort.bam> 

java -jar picard.jar MarkDuplicates
    [-OPTIONS]
    I=<>
    O=<>
    M=<>

java -jar picard.jar MarkDuplicates
    I=<>
    O=<>
    M=<>
   
python cluster_rmdup.py \
    [-OPTIONS]
    <> \
    <>
```

By now you will have an analysis-ready file, free of duplicates.