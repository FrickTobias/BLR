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
test $m == 7f0d3e8bd8d5f147c1f4c4234b4e1efe
