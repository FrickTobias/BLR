#!/bin/bash
set -xeuo pipefail

samtools --version
bowtie2 --version
minimap2 --version
cutadapt --version
starcode --version
snakemake --version
blr --version
samblaster --version
sambamba --version
ema

( cd testdata && bwa index chr1mini.fasta )
( cd testdata && bowtie2-build chr1mini.fasta chr1mini.fasta > /dev/null )
( cd testdata && samtools faidx chr1mini.fasta )

pytest -v tests/

rm -rf outdir-bowtie2
blr init --r1=testdata/reads.1.fastq.gz outdir-bowtie2
blr config \
    --file outdir-bowtie2/blr.yaml \
    --set genome_reference ../testdata/chr1mini.fasta \
    --set reference_variants ../testdata/HG002_GRCh38_GIAB_highconf.chr1mini.vcf \
    --set phasing_ground_truth ../testdata/HG002_GRCh38_GIAB_highconf_triophased.chr1mini.vcf \
    --set max_molecules_per_bc 1

pushd outdir-bowtie2
blr run
m=$(samtools view mapped.sorted.tag.mkdup.bcmerge.mol.filt.bam | md5sum | cut -f1 -d" ")
test $m == 92f61fee39508beadb96f018a6ceac49

# Test phasing
blr run phasing_stats.txt

# Cut away columns 2 and 3 as these change order between linux and osx
m2=$(cut -f1,4- mapped.sorted.tag.mkdup.bcmerge.mol.filt.phase | md5sum | cut -f1 -d" ")
test $m2 == 70c907df8a996d2b3ba3f06fb942b244
