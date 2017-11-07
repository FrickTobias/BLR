# WGH_Analysis

This is a pipeline to perform scaffolding on WGH data. It uses numerous bioinformatics tools with some custom make
parsing tools for formatting in between steps. 

TL;DR: Look at flowchart at the bottom of the page to see pipeline.

## Prerequisites

To run the pipeline, the following software need to be installed:

  - cd-hit-454
  - cutadapt
  - UMItools
  - bowtie2
  
This can be done by writing the following command in your terminal.

```
sudo bash prerequisites.sh
```

It will also be required to have downloaded the programs below. 

  - [fragScaff](https://sourceforge.net/projects/fragscaff/files/?source=navbar) (download fragScaff.pl)
  - [Picard Tools](https://github.com/broadinstitute/picard)

Lastly, a Bowtie2 reference genome (e.g. GR38) is needed, available at Illuminas iGenomes.

[Illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html)

## Useage

First, download the github repository by writing the cloning command in your terminal.

```
git clone https://github.com/FrickTobias/WGH_Analysis.git
```

### Automated analysis
The whole pipeline can be run be using the automation script instead of running all commands individually using 
WGH_Analysis.sh. First however the script need to know where non-pathed software/folders is located on your computer, (Picard 
Tools, bowtie2 reference, fragScaff and the WGH folder) which can be set using setpath.sh.

```
bash setpath.sh <picard_path> <bowtie2_reference> <fragScaff_path>
```
Then the analysis can be started with the following command:
```
bash WGH_automation.sh <read_1.fq> <read_2.fq> <output>
```

### Step by step
Here, every step will be explained in detail. For more details about in-house option useage, see example folder.

#### Read trimming and barcode identification
This is comprised of three steps, all working on read 1. First handle 1 is trimmed from the 5' end and reads with 
those handles are kept. Then the next 20 bases are moved from the sequence up to the header using UMItools extract 
followed by another trimming of handle 2 as in the first step (this only need to be run once to generate paths.sh).

```
cutadapt [-OPTIONS] <read_1.fq> <> <>
UMItools [-OPTIONS] <read_1.fq> <read_2.fq>
cutadapt [-OPTIONS] <read_1.fq> <>
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
python cdhit_prep [-OPTIONS] <read_1.fq> <output_barcode_directory>
for index_file in <barcode_directory>/*; do;
    cd-hit-454 [-OPTIONS] -i $index_file -o $index_file'.clustered'; 
    done
cat <barcode>/.clstr > barcodes.clstr
```

#### Mapping

Mapping is done with bowtie2.

```
bowtie2 [-OPTIONS] <read_1.fq> <read_2.fq> <reference/Bowtie2/genome> <output>
```

#### Scaffolding
The first step is to tag the bam file with the clustering information by using tag_bam.py. After this, picardtools
can be used to remove duplicates (whilst taking barcodes into account). By now, fragScaff can be run twice, first
to estimate N90 of the data and then to perform the final scaffolding.

```
python tag_bam.py [-OPTIONS] <mapped_inserts.bam> <barcodes.clstr> <output_file>
picardtools rmdup
perl fragScaff [-OPTIONS] <mapped_inserts_tagged.bam> <output>
```

## Overview

![](https://github.com/FrickTobias/WGH_Analysis/blob/master/figures/flowchart.png)

_**Figure 1:** Overview. a) Reads are trimmed for sequence NNNNNNNNN, the first handle followed by the following 20 
bases extracted to the header using UMItools extract. Then the second handle is trimmed as the first (but looking 
for the sequence NNNNNNNN). b1) cdhit_prep.py takes the r1_insert.fq, and writes the barcodes to individual files 
according to their indexing bases. b2) CD-HIT-454 is used to cluster all .fa files in the output directory from 
cdhit_prep.py. b3) All .clstr files are concatenated to a NN.clstr file. c) r1 and r2 insert files are mapped using 
bowtie2. d) Mapping files are tagged with cluster number according to NN.clstr file using tag_bam.py in their ‘RG’ tag. 
e) picardtools are used to remove duplicates (taking RG into account.). f) fragScaff is run twice, once to calculate 
N90 of input scaffolds and one more time to create final scaffolds._
