"""
Rules to process 10xGenomics Chromium linked-read raw FASTQs.

READ1 LAYOUT

5'-NNNNNNNNNNNNNNNN NNNNNNN NNN...NNN-3'
   <----barcode---> <-UMI-> <-gDNA-->
         16 nt        7 nt

READ2 LAYOUT

5'-NNNN...NNNN-3'
   <---gDNA-->

Processing is partly based on the 10x end-to-end workflow described for the EMA aligner. See the docs on their github
https://github.com/arshajii/ema#end-to-end-workflow-10x
"""


rule link_to_whitelist:
    output: "barcodes_whitelist.txt"
    shell: "ln -s {config[barcode_whitelist]} {output}"

rule count_10x:
    "Create list of per-barcode count"
    output:
        counts_ncnt = temp("reads.ema-ncnt"),
        counts_fcnt = temp("reads.ema-fcnt")
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
    "Trim reads and bin reads containing the same barcode together. Reads missing barcodes outputed to ema-nobc."
    output:
        bins = directory("temp_bins")
    input:
        r1_fastq="reads.1.fastq.gz",
        r2_fastq="reads.2.fastq.gz",
        counts_ncnt = "reads.ema-ncnt",
        counts_fcnt = "reads.ema-fcnt",
        whitelist = "barcodes_whitelist.txt"
    log: "ema_preproc.log"
    threads: 20
    run:
        hamming_correction = "" if not config["apply_hamming_correction"] else " -h"
        shell(
            "paste <(pigz -c -d {input.r1_fastq} | paste - - - -) <(pigz -c -d {input.r2_fastq} | paste - - - -) |"
            " tr '\t' '\n' |"
            " ema preproc"
            " -w {input.whitelist}"
            " -n {threads}"
            "{hamming_correction}"
            " -t {threads}"
            " -b"
            " -h"
            " -o {output.bins} {input.counts_ncnt} 2>&1 | tee {log}"
        )

rule merge_bins_and_split_pairs:
    "Merge bins of trimmed and barcoded reads together and split into read pairs."
    output:
        r1_fastq="trimmed.barcoded.1.fastq.gz",
        r2_fastq="trimmed.barcoded.2.fastq.gz"
    input:
        bins = "temp_bins"
    shell:
        "cat {input.bins}/ema-bin* |"
        " paste - - - - - - - - |"
        " tee >(cut -f 1-4 | tr '\t' '\n' | pigz -c > {output.r1_fastq}) |"
        " cut -f 5-8 | tr '\t' '\n' | pigz -c > {output.r2_fastq}"

rule split_nobc_reads:
    "Split non-barcoded reads into read pairs."
    output:
        r1_fastq="trimmed.non_barcoded.1.fastq.gz",
        r2_fastq="trimmed.non_barcoded.2.fastq.gz",
    input:
        fastq = "temp_bins/ema-nobc"
    shell:
        " paste - - - - - - - - < {input} |"
        " tee >(cut -f 1-4 | tr '\t' '\n' | pigz -c > {output.r1_fastq}) |"
        " cut -f 5-8 | tr '\t' '\n' | pigz -c > {output.r2_fastq}"
