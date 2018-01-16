# Semi automated analysis
This document describes how to run a semi automated analysis. For examples, see below
 
## Description 

### Read trimming
 
```
bash WGH_read_processing.sh <options> <read1.fq> <read2.fq> <output_dir>
```

### Read mapping

```
bash WGH_mapper.sh <options> <read1.trimmed.fq> <read2.trimmed.fq> <output_bam>
```

### Barcode clustering

```
bash WGH_cluster_bc.sh <options> <read1.trimmed.fq> <bamfile_to_tag> <cluster_output_dir>
```

### Read filtering

```
bash WGH_filtering.sh <options> <mapped.bam> <final_output_file.bam>
```

## Examples

Below an example run is described.

```
bash WGH_read_processing.sh -p 10 -m tobias.frick@scilifelab.se rawdata/180114_dataset_A/A_R1.fq rawdata/180114_dataset_A/A_R2.fq analysis/180114_A;

bash WGH_mapper.sh -p 10 -m tobias.frick@scilifelab.se analysis/180114_A/A_R1.trimmed.fq.gz analysis/180114_A/A_R2.trimmed.fq analysis/180114_A/A.bam

bash WGH_cluster_bc.sh -p 10 -m tobias.frick@scilifelab.se analysis/180114_A/A_R1.trimmed.fq.gz analysis/180114_A/A.bam analysis/180114_A/barcode_clusters 

bash WGH_filtering.sh -p 10 -m tobias.frick@scilifelab.se analysis/180114_A/A.bam analysis/180114_A/A.filt.bam
```