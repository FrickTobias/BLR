#!/bin/bash
set -xeuo pipefail

( cd testdata && bowtie2-build chr1mini.fasta chr1mini > /dev/null )

bash BLR_automation.sh testdata/reads.1.fastq.gz testdata/reads.2.fastq.gz outdir
