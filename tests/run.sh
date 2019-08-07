#!/bin/bash
set -xeuo pipefail

( cd testdata && bowtie2-build chr1mini.fasta chr1mini > /dev/null )

# test single core execution
rm -rf outdir
mkdir -p outdir
ln -s $PWD/testdata/reads.1.fastq.gz outdir/reads.1.fastq.gz
ln -s $PWD/testdata/reads.2.fastq.gz outdir/reads.2.fastq.gz

snakemake --configfile tests/test_config.yaml outdir/reads.1.final.fastq.gz \
    outdir/reads.2.final.fastq.gz

m=$(samtools view outdir/mapped.sorted.tag.rmdup.x2.filt.bam | md5sum | cut -f1 -d" ")
test $m == 424db26031fa6d874da95ad5f023bbce

# test parallel execution
rm -rf outdir_para
mkdir -p outdir_para
ln -s $PWD/testdata/reads.1.fastq.gz outdir_para/reads.1.fastq.gz
ln -s $PWD/testdata/reads.2.fastq.gz outdir_para/reads.2.fastq.gz

snakemake --configfile tests/test_config.yaml outdir_para/reads.1.final.fastq.gz \
    outdir_para/reads.2.final.fastq.gz -j 3

# Check depth coverage instead of sam file as the BC index will change in parallel execution.
m2=$(samtools depth outdir_para/mapped.sorted.tag.rmdup.x2.filt.bam | md5sum | cut -f1 -d" ")
test $m2 == 8fcc72b9d2d6ff8510e55d5d36937626

# Check that depth is the same for the single core file.
m3=$(samtools depth outdir/mapped.sorted.tag.rmdup.x2.filt.bam | md5sum | cut -f1 -d" ")
test $m2 == $m3