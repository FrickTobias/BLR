# kate: syntax Python;
rule:
    "Trim away E handle on R1 5'. Also removes reads shorter than 85 bp."
    output:
        r1_fastq="{dir}/trimmed-a.1.fastq.gz",
        r2_fastq="{dir}/trimmed-a.2.fastq.gz"
    input:
        r1_fastq="{dir}/reads.1.fastq.gz",
        r2_fastq="{dir}/reads.2.fastq.gz"
    log: "{dir}/trimmed-a.log"
    threads: 20
    shell:
        "cutadapt -g ^CAGTTGATCATCAGCAGGTAATCTGG"
        " -j {threads}"
        " --discard-untrimmed"
        " -e 0.2"
        " -m 65"
        " -o {output.r1_fastq}"
        " -p {output.r2_fastq}"
        " {input.r1_fastq}"
        " {input.r2_fastq}"
        " | tee {log}"
