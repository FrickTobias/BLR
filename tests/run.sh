#!/bin/bash
set -xeuo pipefail

( cd testdata && bowtie2-build chr1mini.fasta chr1mini > /dev/null )

rm -rf outdir
mkdir -p outdir
ln -s $PWD/testdata/reads.1.fastq.gz outdir/reads.1.fastq.gz
ln -s $PWD/testdata/reads.2.fastq.gz outdir/reads.2.fastq.gz

snakemake --configfile tests/test_config.yaml outdir/reads.1.final.fastq.gz \
    outdir/reads.2.final.fastq.gz

m=$(samtools sort -n outdir/mapped.sorted.tag.mkdup.bcmerge.filt.bam | samtools view - | md5sum | cut -f1 -d" ")
test $m == be2dbe2f5a2ab660a949e1944e077a79

# Test phasing
snakemake --configfile tests/test_config.yaml outdir/mapped.sorted.tag.mkdup.bcmerge.filt.phase \
    outdir/mapped.sorted.tag.mkdup.bcmerge.filt.phase.phased.vcf

m2=$(md5sum outdir/mapped.sorted.tag.mkdup.bcmerge.filt.phase | cut -f1 -d" ")
test $m2 == 5959009bb88bee382932a4080954cdb3
