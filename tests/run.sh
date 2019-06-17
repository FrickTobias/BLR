#!/bin/bash
set -xeuo pipefail

( cd testdata && bowtie2-build chr1mini.fasta chr1mini > /dev/null )

rm -rf outdir
bash BLR_automation.sh testdata/reads.1.fastq.gz testdata/reads.2.fastq.gz outdir
m=$(samtools view outdir/reads.1.fastq.sort.tag.rmdup.x2.filt.bam | md5 | cut -f1 -d" ")
test $m == 51788863e44c1a43f8da7287b1cb2867
