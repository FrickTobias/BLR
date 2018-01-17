# WGH_Analysis

This is a pipeline for handling of [WGH]() data. It takes fastq files as input and outputs a bam file ready for variant calling. For a quick overview, look at the [flowchart](https://github.com/FrickTobias/WGH_Analysis/blob/master/README.md#overview) to see what the pipeline does.

## Prerequisites

To run the pipeline, the following software need to be installed:

  - [cd-hit-454](https://github.com/weizhongli/cdhit.git)
  - [cutadapt](https://github.com/marcelm/cutadapt.git)
  - [bowtie2](https://github.com/BenLangmead/bowtie2)
  - [samtools](https://github.com/samtools/samtools)
  - [qvalues](https://github.com/jeffhsu3/qvalue.git)
  
This can be done by writing the following command in your terminal which will install everything but qvalues (since pip 
will default to the python 2 version). 

```
sudo bash prerequisites.sh
```

To get qvalues for python 3, clone the github fork above and install it.

```
git clone https://github.com/jeffhsu3/qvalue.git
```
```
python3 setup.py build
python3 setup.py install
```

It will also be required to have downloaded [Picard Tools](https://github.com/broadinstitute/picard) and a Bowtie2 reference genome (e.g. GRCh38), available at e.g. [Illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html).

## Useage

### Setup

First, download the github repository by writing the cloning command in your terminal.

```
git clone https://github.com/FrickTobias/WGH_Analysis.git
```

Then provide WGH_Analysis with the appropriate paths for Picard Tools and your Bowtie2 reference data. 


```
bash setpath.sh <picard_path> <bowtie2_reference>
```

### Automated analysis
The whole pipeline can be run be using the automation script WGH_Analysis.sh. For standard useage, run the following command.

```
bash WGH_automation.sh -m <john.doe@myworkplace.com -p <processors> <read_1.fq> <read_2.fq> <output> 
```

For all available options, see -h (--help) and for more details consult the [step-by-step](https://github.com/FrickTobias/WGH_Analysis/tree/master/step-by-step) 
folder which describes all steps performed by WGH_automation. For examples, see the [example folder](https://github.com/FrickTobias/WGH_Analysis/tree/master/example) 
where an example run is thouroughly described.

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
