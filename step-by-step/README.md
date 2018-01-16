### Main steps
Here, every step will be explained in detail. For more details about in-house option useage, see example folder. There
are also semi-automated scripts for performing the main parts of the analysis, see semi_auto.md.

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

#### Filtering

