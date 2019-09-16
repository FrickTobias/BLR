#!/bin/bash
set -xeuo pipefail

( cd testdata && bowtie2-build chr1mini.fasta chr1mini > /dev/null )

rm -rf outdir
mkdir -p outdir
ln -s $PWD/testdata/reads.1.fastq.gz outdir/reads.1.fastq.gz
ln -s $PWD/testdata/reads.2.fastq.gz outdir/reads.2.fastq.gz
cp tests/test_config.yaml outdir/config.yaml
cd outdir
blr run reads.1.final.fastq.gz reads.2.final.fastq.gz

m=$(samtools sort -n mapped.sorted.tag.mkdup.bcmerge.mol.filt.bam | samtools view - | md5sum | cut -f1 -d" ")
test $m == 134db15680443fc32d25ba577a0c5700

# Test phasing
blr run phasing_stats.txt

# Cut away columns 2 and 3 as these change order between linux and osx
m2=$(cut -f1,4- mapped.sorted.tag.mkdup.bcmerge.mol.filt.phase | md5sum | cut -f1 -d" ")
test $m2 == 70c907df8a996d2b3ba3f06fb942b244
