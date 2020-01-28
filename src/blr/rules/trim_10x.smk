"""
Rules to process 10xGenomics Chromium linked-read raw FASTQs.

READ1 LAYOUT

5'-NNNNNNNNNNNNNNNN NNNNNNN NNN...NNN-3'
   <----barcode---> <-UMI-> <-gDNA-->
         16 nt        7 nt

READ2 LAYOUT

5'-NNNN...NNNN-3'
   <---gDNA-->

"""


rule link_to_whitelist:
    output: "barcodes_whitelist.txt"
    shell: "ln -s {config[barcode_whitelist]} {output}"

rule count_10x:
    output:
        counts_ncnt = temp("reads.ema-ncnt"),
        counts_fcnt = temp("reads.ema-fcnt"),
    input:
        r1_fastq="reads.1.fastq.gz",
        r2_fastq="reads.2.fastq.gz",
        whitelist="barcodes_whitelist.txt"
    log: "ema_count.log"
    shell:
        "paste <(pigz -c -d {input.r1_fastq} | paste - - - -) <(pigz -c -d {input.r2_fastq} | paste - - - -) |"
        " tr '\t' '\n' |"
        " ema count"
        " -w {input.whitelist}"
        " -o reads 2> {log}"

rule preproc_10x:
    output:
        bins = temp(directory("temp_bins"))
    input:
        r1_fastq="reads.1.fastq.gz",
        r2_fastq="reads.2.fastq.gz",
        counts_ncnt = "reads.ema-ncnt",
        counts_fcnt = "reads.ema-fcnt",
        whitelist = "barcodes_whitelist.txt"
    log: "ema_preproc.log"
    threads: 20
    shell:
        "paste <(pigz -c -d {input.r1_fastq} | paste - - - -) <(pigz -c -d {input.r2_fastq} | paste - - - -) |"
        " tr '\t' '\n' |"
        " ema preproc"
        " -w {input.whitelist}"
        " -n {threads}"
        " -t {threads}"
        " -b"
        " -h"
        " -o {output.bins} {input.counts_ncnt} 2>&1 | tee {log}"

rule merge_bins_and_split_pairs:
    output:
        r1_fastq="trimmed.barcoded.1.fastq.gz",
        r2_fastq="trimmed.barcoded.2.fastq.gz",

    input:
        bins = "temp_bins"
    shell:
        "cat {input.bins}/ema-bin* |"
        " paste - - - - - - - - |"
        " tee >(cut -f 1-4 | tr '\t' '\n' | pigz -c > {output.r1_fastq}) |"
        " cut -f 5-8 | tr '\t' '\n' | pigz -c > {output.r2_fastq}"

