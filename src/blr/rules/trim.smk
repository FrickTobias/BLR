"""
Rules connected to trimming of fastq files.
"""

rule trim:
    # Trim away E handle on R1 5'. Also removes reads shorter than 85 bp.
    # Extract barcode sequence and place in header.
    # Cut H1691' + TES sequence from 5' of R1. H1691'=CATGACCTCTTGGAACTGTC, TES=AGATGTGTATAAGAGACAG.
    # Cut 3' TES' sequence from R1 and R2. TES'=CTGTCTCTTATACACATCT
    # Discard untrimmed.
    output:
        r1_fastq="trimmed-c.1.fastq.gz",
        r2_fastq="trimmed-c.2.fastq.gz"
    input:
        r1_fastq="reads.1.fastq.gz",
        r2_fastq="reads.2.fastq.gz"
    log:
        a = "trimmed-a.log",
        b = "trimmed-b.log"
    threads: 20
    shell:
        "cutadapt" #Initial trim
        " -g 'XNNNCAGTTGATCATCAGCAGGTAATCTGG;min_overlap=26'"
        " -e 0.2"
        " --discard-untrimmed"
        " -j {threads}"
        " -m 65"
        " -o -"
        " --interleaved"
        " {input.r1_fastq}"
        " {input.r2_fastq}"
        " 2> {log.a} |"
        " blr extractbarcode" #Extract barcodes
        " -"
        " -o1 - |" 
        " cutadapt" # Final trim
        " -a ^CATGACCTCTTGGAACTGTCAGATGTGTATAAGAGACAG...CTGTCTCTTATACACATCT"
        " -A CTGTCTCTTATACACATCT"
        " -e 0.2"
        " --discard-untrimmed"
        " --pair-filter 'first'"
        " -j {threads}"
        " -m 25"
        " -o {output.r1_fastq}"
        " -p {output.r2_fastq}"
        " -"
        " --interleaved"
        " > {log.b}"
