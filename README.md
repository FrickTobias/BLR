# WGH_Analysis

Description
  - softwares
  - custom-written parsing scripts

TL;DR: Look at overview.

## Prerequisites

To install the required software for this pipeline, write the following in your terminal.
```
sudo apt-get install prerequisites.txt
```
This will install (fragScaff, cdhit_prep and tag_bam.py scripts do not require installation and can be run by writing perl/python before the
softwares name)
  - cd-hit-45
  - cutadapt
  - UMItools
  - bowtie2
  
Furthermore you will also need a Bowtie2 reference genome, which can be downloaded from NCBI (link below for GR38)
```
Bowtie2 reference (e.g. GR38)
```
## Useage

### Automated analysis
WGH_automation description

### Step by step
Every step with description of options

#### Read trimming & barcode identification
description

```
cutadapt first_handle keep_reads
UMItools extract bc
cutadapt second_handle keep_reads
```

#### Clustering
description
```
python cdhit_prep -r 2 -f 0
cd-hit-454 AT.fa OPTIONS
cat folder/.clstr > NN.clstr
```

#### Mapping
Description
```
bowtie2 insert reference/Bowtie2/genome
```

#### Scaffolding
Description - mainly about why need to tag bam file.
```
python tag_bam.py bam clstr output
fragscaff
```

## Overview
See project wiki, link below.

![](https://github.com/FrickTobias/WGH_Analysis/blob/master/figures/flowchart.png)

Figure 1: Overview. a) Reads are trimmed for sequence NNNNNNNNN, the first handle followed by the following 20 bases 
extracted to the header using UMItools extract. Then the second handle is trimmed as the first (but looking for the 
sequence NNNNNNNN). b1) cdhit_prep.py takes the r1_insert.fq, and writes the barcodes to individual files according 
to their indexing bases. b2) CD-HIT-454 is used to cluster all .fa files in the output directory from cdhit_prep.py. 
b3) All .clstr files are concatenated to a NN.clstr file. c) r1 and r2 insert files are mapped using bowtie2. d) 
Mapping files are tagged with cluster number according to NN.clstr file using tag_bam.py in their ‘RG’ tag. e) 
picardtools are used to remove duplicates (taking RG into account.). f) fragScaff is run twice, once to calculate 
N90 of input scaffolds and one more time to create final scaffolds.

