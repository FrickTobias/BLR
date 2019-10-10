#!/bin/bash
set -xeuo pipefail

samtools --version
bowtie2 --version
picard SamToFastq --version || true
cutadapt --version
starcode --version
snakemake --version
blr --version

( cd testdata && bwa index chr1mini.fasta )

rm -rf outdir-bwa
blr init --r1=testdata/reads.1.fastq.gz outdir-bwa
sed 's|read_mapper: .*|read_mapper: bwa|' tests/test_config.yaml > outdir-bwa/blr.yaml

pushd outdir-bwa
blr run
m=$(samtools view mapped.sorted.tag.mkdup.bcmerge.mol.filt.bam | md5sum | cut -f1 -d" ")
test $m == b7bffb901d59030ea4cc939a29c85643
popd

( cd testdata && bowtie2-build chr1mini.fasta chr1mini.fasta > /dev/null )

rm -rf outdir-bowtie2
blr init --r1=testdata/reads.1.fastq.gz outdir-bowtie2
cp tests/test_config.yaml outdir-bowtie2/blr.yaml
pushd outdir-bowtie2
blr run
m=$(samtools view mapped.sorted.tag.mkdup.bcmerge.mol.filt.bam | md5sum | cut -f1 -d" ")
test $m == 5607178a324ce4a394e5370a4d192377

# Test phasing
blr run phasing_stats.txt

# Cut away columns 2 and 3 as these change order between linux and osx
m2=$(cut -f1,4- mapped.sorted.tag.mkdup.bcmerge.mol.filt.phase | md5sum | cut -f1 -d" ")
test $m2 == 70c907df8a996d2b3ba3f06fb942b244
