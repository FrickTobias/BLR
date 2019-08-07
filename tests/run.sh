#!/bin/bash
set -xeuo pipefail

( cd testdata && bowtie2-build chr1mini.fasta chr1mini > /dev/null )

# test single core execution
rm -rf outdir
mkdir -p outdir
ln -s $PWD/testdata/reads.1.fastq.gz outdir/reads.1.fastq.gz
ln -s $PWD/testdata/reads.2.fastq.gz outdir/reads.2.fastq.gz

snakemake --configfile tests/test_config.yaml outdir/reads.1.final.fastq.gz \
    outdir/reads.2.final.fastq.gz -j 3

m=$(samtools sort -n outdir/mapped.sorted.tag.rmdup.x2.filt.bam | samtools view - | md5sum | cut -f1 -d" ")
test $m == ed01044fce9984087c5a8f015e9e1bd2