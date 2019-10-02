"""
Rules connected to trimming of FASTQ files.

READ1 LAYOUT

5'-CAGTTGATCATCAGCAGGTAATCTGGBDVHBDVHBDVHBDVHBDVHCATGACCTCTTGGAACTGTCAGATGTGTATAAGAGACAGNNNN...NNNN(CTGTCTCTTATACACATCT)-3'
   <------------h1----------><-------DBS--------><-----------------h2------------------><---gDNA--><---------h3-------->
    h1 may inlude frameshift                                                                        Presence depencd on
    oligos varying from 0-4                                                                         insert length
    extra oligos in 5' end

READ 2 LAYOUT

5'-NNNN...NNNN(CTGTCTCTTATACACATCT)-3'
   <---gDNA--><---------h3-------->
              Presence depencd on
              insert length
"""

h1="CAGTTGATCATCAGCAGGTAATCTGG"
DBS="N"*20
h2="CATGACCTCTTGGAACTGTCAGATGTGTATAAGAGACAG"
h3="CTGTCTCTTATACACATCT"
trim_len=sum(map(len, [h1,DBS,h2]))
extract_len=len(h1)

rule trim:
    # Trim away 5' and possible 3' handles on read1 and trim possible 3' handles on read2.
    output:
        r1_fastq="trimmed.1.fastq.gz",
        r2_fastq="trimmed.2.fastq.gz"
    input:
        r1_fastq="reads.1.fastq.gz",
        r2_fastq="reads.2.fastq.gz",
    log: "cutadapt_trim.log"
    threads: 20
    shell:
        "cutadapt" #Initial trim
        " -g 'XNNN{h1}{DBS}{h2};min_overlap={trim_len}...{h3};optional'"
        " -A {h3}"
        " --pair-filter 'first'"
        " -e 0.2"
        " --discard-untrimmed"
        " -j {threads}"
        " -m 25"
        " -o {output.r1_fastq}"
        " -p {output.r2_fastq}"
        " {input.r1_fastq}"
        " {input.r2_fastq}"
        "> {log} "

rule extract_DBS:
    # Extract barcode sequence from read1 FASTQ
    output:
        fastq="barcodes.fastq",
    input:
        fastq="reads.1.fastq.gz",
    log: "cutadapt_extract_DBS.log"
    threads: 20
    shell:
        "cutadapt" #Initial trim
        " -g 'XNNN{h1};min_overlap={extract_len}...{h2}'"
        " -e 0.2"
        " --discard-untrimmed"
        " -j {threads}"
        " -m 19"
        " -M 21"
        " -o {output.fastq}"
        " {input.fastq}"
        " > {log}"
