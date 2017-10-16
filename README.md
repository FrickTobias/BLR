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
```
[WGH Analysis flowchart](https://github.com/FrickTobias/WGH_Analysis/wiki/Flowchart)
```

