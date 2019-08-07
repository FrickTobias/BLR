from snakemake.utils import available_cpu_count
#
# Piped trimming of reads, executed if -j/--jobs/--cores is given with number larger
# then or equal to 3 as there are 3 processes in the pipe. Intermediate fastqs are interleaved.
#
rule trim_r1_handle_pipe:
    #Trim away E handle on R1 5'. Also removes reads shorter than 85 bp.
    output:
        interleaved_fastq=pipe("{dir}/trimmed-a.fastq")
    input:
        r1_fastq="{dir}/reads.1.fastq.gz",
        r2_fastq="{dir}/reads.2.fastq.gz"
    log: "{dir}/trimmed-a.log"
    threads: int(available_cpu_count()/3)
    shell:
        "cutadapt"
        " -g ^CAGTTGATCATCAGCAGGTAATCTGG"
        " -e 0.2"
        " --discard-untrimmed"
        " -j {threads}"
        " -m 65"
        " -o {output.interleaved_fastq}"
        " --interleaved "
        " {input.r1_fastq}"
        " {input.r2_fastq}"
        " > {log}"

rule extract_barcodes_pipe:
    output:
        interleaved_fastq=pipe("{dir}/unbarcoded.fastq")
    input:
        interleaved_fastq="{dir}/trimmed-a.fastq"
    shell:
        # BDHVBDVHBDVHBDVH
        "blr extractbarcode"
        " {input.interleaved_fastq}"
        " -o1 {output.interleaved_fastq}"

rule final_trim_pipe:
    # Cut H1691' + TES sequence from 5' of R1. H1691'=CATGACCTCTTGGAACTGTC, TES=AGATGTGTATAAGAGACAG.
    # Cut 3' TES' sequence from R1 and R2. TES'=CTGTCTCTTATACACATCT
    # Discard untrimmed.
    output:
        r1_fastq="{dir}/trimmed-c.1.fastq.gz",
        r2_fastq="{dir}/trimmed-c.2.fastq.gz"
    input:
        interleaved_fastq="{dir}/unbarcoded.fastq"
    log: "{dir}/trimmed-b.log"
    threads: int(available_cpu_count()/3)
    shell:
        "cutadapt"
        " -a ^CATGACCTCTTGGAACTGTCAGATGTGTATAAGAGACAG...CTGTCTCTTATACACATCT "
        " -A CTGTCTCTTATACACATCT "
        " -e 0.2"
        " --discard-untrimmed"
        " --pair-filter 'first'"
        " -j {threads}"
        " -m 25"
        " -o {output.r1_fastq}"
        " -p {output.r2_fastq}"
        " {input.interleaved_fastq}"
        " --interleaved"
        " > {log}"
