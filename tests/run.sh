#!/bin/bash
set -xeuo pipefail

( cd testdata && bowtie2-build chr1mini.fasta chr1mini > /dev/null )

rm -rf outdir
bash BLR_automation.sh testdata/reads.1.fastq.gz testdata/reads.2.fastq.gz outdir
m=$(samtools view outdir/reads.1.fastq.sort.tag.rmdup.x2.filt.bam | md5sum | cut -f1 -d" ")
test $m == 424db26031fa6d874da95ad5f023bbce