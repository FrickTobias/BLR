# Barcode_activity

Description
  - softwares
  - custom-written parsing scripts

## Prerequisites

```
Software list + links
```
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

## Examples
### Automation analalysis
```
bash WGH_automation ...
```

### Step by step
```
cutadapt ...
UMItools extract
cutadadpt
```

```
cdhit_prep.py
cdhit
cat * > NN.fq
```
```
...
```

## Overview

Flowchart
(toned box r1/r2.fq => small .filetype => boxed software funtions => ...)
