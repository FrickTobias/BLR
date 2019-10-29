#!/bin/bash
set -xeuo pipefail

samtools --version
bowtie2 --version
minimap2 --version
picard SamToFastq --version || true
cutadapt --version
starcode --version
snakemake --version
blr --version
samblaster --version
sambamba --version

( cd testdata && bwa index chr1mini.fasta )
( cd testdata && bowtie2-build chr1mini.fasta chr1mini.fasta > /dev/null )

pytest -v tests/

rm -rf outdir-bowtie2
blr init --r1=testdata/reads.1.fastq.gz outdir-bowtie2
cp tests/test_config.yaml outdir-bowtie2/blr.yaml
pushd outdir-bowtie2
blr run
m=$(samtools view mapped.sorted.tag.mkdup.bcmerge.mol.filt.bam | md5sum | cut -f1 -d" ")
test $m == e389e3f6e8ba14466f232b19fa1a0ee5

# Test phasing
blr run phasing_stats.txt

# Cut away columns 2 and 3 as these change order between linux and osx
m2=$(cut -f1,4- mapped.sorted.tag.mkdup.bcmerge.mol.filt.phase | md5sum | cut -f1 -d" ")
test $m2 == 70c907df8a996d2b3ba3f06fb942b244
