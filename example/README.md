### Examples of how run pipeline

#### Basics

List information about script and options available (plus default settings):

```
bash WGH_automation.sh -h
```

Most basic useage:

```
bash WGH_automation.sh testdata_read1.fastq.gz testdata_read2.fastq.gz testdata_analysis
```

Recommended options, for running on a computer with (-p) 24 processors, (-m) mailing 
to john.doe@domain.com and (-r) removing redundant files created in the pipeline:

#### Advanced useage

```
bash WGH_automation.sh -r -m john.doe@domain.com -p 24 testdata_read1.fastq.gz testdata_read2.fastq.gz testdata_analysis
```

(-s) Starting at step 2 and (-e) ending after step 3 (see -h, pipeline outline for information about more information 
about step numbers):

```
bash WGH_automation.sh -s 2 -e 3 testdata_read1.fastq.gz testdata_read2.fastq.gz testdata_analysis
```

Running with 4 barcode index nucleotides:

```
bash WGH_automation.sh -i 4 testdata_read1.fastq.gz testdata_read2.fastq.gz testdata_analysis
```
