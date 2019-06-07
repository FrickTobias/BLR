# kate: syntax Python;
rule:
    "Trim away E handle on R1 5'. Also removes reads shorter than 85 bp."
    output:
        r1_fastq=temp("{dir}/trimmed-a.1.fastq.gz"),
        r2_fastq=temp("{dir}/trimmed-a.2.fastq.gz")
    input:
        r1_fastq="{dir}/reads.1.fastq.gz",
        r2_fastq="{dir}/reads.2.fastq.gz"
    log: "{dir}/trimmed-a.log"
    threads: 20
    shell:
        "cutadapt"
        " -g ^CAGTTGATCATCAGCAGGTAATCTGG"
        " -e 0.2"
        " --discard-untrimmed"
        " -j {threads}"
        " -m 65"
        " -o {output.r1_fastq}"
        " -p {output.r2_fastq}"
        " {input.r1_fastq}"
        " {input.r2_fastq}"
        " | tee {log}"


rule:
    """Cut TES from 5' of R1. TES=AGATGTGTATAAGAGACAG. Discard untrimmed."""
    output:
        r1_fastq=temp("{dir}/trimmed-b.1.fastq.gz"),
        r2_fastq=temp("{dir}/trimmed-b.2.fastq.gz")
    input:
        r1_fastq="{dir}/unbarcoded.1.fastq.gz",
        r2_fastq="{dir}/unbarcoded.2.fastq.gz"
    log: "{dir}/trimmed-b.log"
    threads: 20
    shell:
        "cutadapt"
        " -g AGATGTGTATAAGAGACAG"
        " -e 0.2"
        " --discard-untrimmed"
        " -j {threads}"
        " -o {output.r1_fastq}"
        " -p {output.r2_fastq}"
        " {input.r1_fastq}"
        " {input.r2_fastq}"
        " | tee {log}"
