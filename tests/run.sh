#!/bin/bash
set -xeuo pipefail

samtools --version
bowtie2 --version
picard SamToFastq --version || true
cutadapt --version
cd-hit --help | head -n 1 || true
snakemake --version
blr --version

( cd testdata && bowtie2-build chr1mini.fasta chr1mini > /dev/null )

rm -rf outdir
blr init --r1=testdata/reads.1.fastq.gz outdir
cp tests/test_config.yaml outdir/config.yaml
cd outdir
blr run

m=$(samtools view mapped.sorted.tag.mkdup.bcmerge.mol.filt.bam | md5sum | cut -f1 -d" ")
test $m == 827612d1dd59d07071defa26bc8add4c

# Test phasing
blr run phasing_stats.txt

# Cut away columns 2 and 3 as these change order between linux and osx
m2=$(cut -f1,4- mapped.sorted.tag.mkdup.bcmerge.mol.filt.phase | md5sum | cut -f1 -d" ")
test $m2 == 70c907df8a996d2b3ba3f06fb942b244
