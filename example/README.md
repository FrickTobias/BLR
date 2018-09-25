## Examples of how run pipeline

### Basics

List information about script and options available (plus default settings):

```
bash BLR_automation.sh -h
```

Most basic useage, which runs the pipeline using default options and writes output to a new folder called 
testdata_analysis_folder_output, found under example folder:

```
bash BLR_automation.sh example/testdata_read1.fastq.gz example/testdata_read2.fastq.gz example/testdata_analysis_folder_output
```

Recommended options, for running on a computer with (-p) 24 processors, (-m) mailing 
to john.doe@myworkplace.com when completing pipeline steps and (-r) removing software log files created in the pipeline:

```
bash BLR_automation.sh -r -m john.doe@myworkplace.com -p 24 testdata_read1.fastq.gz testdata_read2.fastq.gz testdata_analysis
```

### Advanced useage

(-s) Starting at step 2 and (-e) ending after step 3 (see -h, pipeline outline for information about more information 
about step numbers):

```
bash BLR_automation.sh -s 2 -e 3 testdata_read1.fastq.gz testdata_read2.fastq.gz testdata_analysis
```

Running with 4 barcode index nucleotides:

```
bash BLR_automation.sh -i 4 testdata_read1.fastq.gz testdata_read2.fastq.gz testdata_analysis
```

### Testdata output

Command:

```
bash BLR_automation.sh -f -r -m john.doe@myworkplace.com -p 24 testdata_read1.fastq.gz testdata_read2.fastq.gz testdata_analysis
```

Output:

```
0. Argparsing & options
Read 1:		testdata_read1.fastq.gz
Read 2:		testdata_read2.fastq.gz
Output:		testdata_analysis
 
Threads:	    24
Starts at step:	1
End at step:	4
Mail:           john.doe@myworkplace.com
 
Fri Mar  2 11:48:35 CET 2018	ANALYSIS STARTING
 
1. Demultiplexing
Fri Mar  2 11:48:35 CET 2018	1st adaptor removal
Fri Mar  2 11:48:36 CET 2018	1st adaptor removal done
Fri Mar  2 11:48:36 CET 2018	Barcode extraction
Fri Mar  2 11:48:36 CET 2018	Barcode extraction done
Fri Mar  2 11:48:36 CET 2018	2nd adaptor removal
Fri Mar  2 11:48:36 CET 2018	2nd adaptor removal done
Fri Mar  2 11:48:36 CET 2018	3' trimming
Fri Mar  2 11:48:37 CET 2018	3' trimming done
Fri Mar  2 11:48:37 CET 2018	Intact reads: 0.973162 %
 
2. Mapping
Fri Mar  2 11:48:37 CET 2018	Mapping
Fri Mar  2 11:48:39 CET 2018	Mapping done
Fri Mar  2 11:48:39 CET 2018	Sorting
Fri Mar  2 11:48:39 CET 2018	Sorting done
Fri Mar  2 11:48:39 CET 2018	1st map stats
Fri Mar  2 11:48:39 CET 2018	1st map stats done
Fri Mar  2 11:48:39 CET 2018	Filtering
Fri Mar  2 11:48:39 CET 2018	Filtering done
Fri Mar  2 11:48:39 CET 2018	2nd map stats
Fri Mar  2 11:48:39 CET 2018	2nd map stats done
 
3. Clustering
Fri Mar  2 11:48:39 CET 2018	Barcode fasta generation
Fri Mar  2 11:48:39 CET 2018	Barcode fasta generation done
Fri Mar  2 11:48:39 CET 2018	Barcode clustering
Fri Mar  2 11:48:40 CET 2018	Barcode clustering done
Fri Mar  2 11:48:40 CET 2018	Bam tagging
Fri Mar  2 11:48:41 CET 2018	Bam tagging done
 
4. Duplicate removal
Fri Mar  2 11:48:41 CET 2018	Duplicate removal
Fri Mar  2 11:48:57 CET 2018	Duplicate removal done
Fri Mar  2 11:48:57 CET 2018	Barcode duplicate marking
Fri Mar  2 11:49:15 CET 2018	Barcode duplicate marking done
Fri Mar  2 11:49:15 CET 2018	Cluster merging
Fri Mar  2 11:49:15 CET 2018	Cluster merging done
Fri Mar  2 11:49:15 CET 2018	Fastq generation
Fri Mar  2 11:49:16 CET 2018	Fastq generation done
 
Fri Mar  2 11:49:16 CET 2018	ANALYSIS FINISHED
```

Generated folder content:

```
1_trim.log						
2_map.log						
3_cluster.log					
4_rmdup.log								
testdata_read1.fastq.trimmed.fastq.gz
testdata_read2.fastq.trimmed.fastq.gz
testdata_read1.fastq.sort.bam
NNN.clstr
testdata_read1.fastq.sort.filt.tag.bam
testdata_read1.fastq.sort.filt.tag.rmdup.mkdup.bam
testdata_read1.fastq.sort.filt.tag.rmdup.x2.bam
testdata_read1.fastq.sort.filt.tag.rmdup.x2.bam.log
testdata_read1.fastq.final.fastq.gz	
testdata_read2.fastq.final.fastq.gz

```

### File contents

Logfiles

```
<FILE>                              <CONTENT>
1_trim.log                          Handle integrity
2_map.log                           Mapping details 
3_cluster.log                       Counts (number unique sequences) for barcode index files
4_rmdup.log                         Counts for duplicate marked reads

```
.fastq files
```					
<FILE>                              <CONTENT>		
R1.trimmed.fastq.gz                 Demultiplexed & trimmed reads
R2.trimmed.fastq.gz                 Demultiplexed & trimmed reads
R1.final.fastq.gz                   Final filtered reads
R2.final.fastq.gz                   Final filtered reads

```
.bam files
```	
<FILE>                              <CONTENT>
.sort.bam                           Mapped inserts
.sort.filt.tag.bam                  -..-, filtered for unmapped reads, non-primary alignments, tagged with barcode sequences
.sort.filt.tag.rmdup.mkdup.bam      -..-, with read duplicates removed and cluster duplicates marked
.sort.filt.tag.rmdup.x2.bam         -..-, with cluster duplicates merged

```
Barcode files
```	
<FILE>                              <CONTENT>
.clstr                              Clustered barcode sequences 
```