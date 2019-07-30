# kate: syntax Python;
import itertools

# Parameters
index_nucleotides = 3
indexes = ["".join(tuple) for tuple in itertools.product("ATCG", repeat=index_nucleotides)] if index_nucleotides > 0 else ["all"]
cluster_tag="BC"


# Currently paths are read from a paths.txt file into a dict, possibly we would want a config file for this.
paths = {}
with open("paths.txt","r") as paths_file:
    for line in paths_file.readlines():
        name, path = line.strip().split("=")
        paths[name] = path

rule trim_r1_handle:
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
        " > {log}"


rule extract_barcodes:
    output:
        r1_fastq=temp("{dir}/unbarcoded.1.fastq"),
        r2_fastq=temp("{dir}/unbarcoded.2.fastq")
    input:
        r1_fastq="{dir}/trimmed-a.1.fastq.gz",
        r2_fastq="{dir}/trimmed-a.2.fastq.gz"
    log: "{dir}/extractbarcode.log"
    shell:
        # BDHVBDVHBDVHBDVH
        "blr extractbarcode"
        " {input.r1_fastq} {input.r2_fastq}"
        " -o1 {output.r1_fastq} -o2 {output.r2_fastq}"
        " 2> {log}"


rule compress:
    output: "{dir}/unbarcoded.{nr}.fastq.gz"
    input: "{dir}/unbarcoded.{nr}.fastq"
    shell:
        "pigz < {input} > {output}"

rule final_trim:
    """
    Cut H1691' + TES sequence from 5' of R1. H1691'=CATGACCTCTTGGAACTGTC, TES=AGATGTGTATAAGAGACAG.
    Cut 3' TES' sequence from R1 and R2. TES'=CTGTCTCTTATACACATCT
    Discard untrimmed.
    """
    output:
        r1_fastq="{dir}/trimmed-c.1.fastq.gz",
        r2_fastq="{dir}/trimmed-c.2.fastq.gz"
    input:
        r1_fastq="{dir}/unbarcoded.1.fastq.gz",
        r2_fastq="{dir}/unbarcoded.2.fastq.gz"
    log: "{dir}/trimmed-b.log"
    threads: 20
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
        " {input.r1_fastq}"
        " {input.r2_fastq}"
        " > {log}"

# TODO add description to new rules.
if index_nucleotides == 0:
    rule cdhitprep:
        output:
            "{dir}/unique_bc/all.fa"
        input:
            r1_fastq = "{dir}/reads.1.fastq.trimmed.fastq.gz"
        log:
            stdout = "{dir}/cdhit_prep.stdout",
            stderr = "{dir}/cdhit_prep.stderr"
        shell:
            "blr cdhitprep "
            " {input.r1_fastq}"
            " {output}"
            " -f 0 > {log.stdout} 2> {log.stderr}"
else:
    rule cdhitprep:
        output:
            expand("{{dir}}/unique_bc/{sample}.fa", sample=indexes)
        input:
            r1_fastq = "{dir}/reads.1.fastq.trimmed.fastq.gz"
        params:
            dir = "{dir}/unique_bc/"
        log:
            stdout = "{dir}/cdhit_prep.stdout",
            stderr = "{dir}/cdhit_prep.stderr"
        shell:
            "blr cdhitprep "
            " {input.r1_fastq}"
            " {params.dir}"
            " -i {index_nucleotides}"
            " -f 0 > {log.stdout} 2> {log.stderr}"

rule barcode_clustering:
    input:
       "{dir}/unique_bc/{sample}.fa"
    output:
        "{dir}/unique_bc/{sample}.clustered",
        "{dir}/unique_bc/{sample}.clustered.clstr"
    threads: 20
    log: "{dir}/unique_bc/{sample}.clustering.log"
    params:
        prefix= lambda wc,output: os.path.splitext(output[1])[0]
    shell:
        " (cd-hit-454 "
        " -i {input} "
        " -o {params.prefix} "
        " -T {threads} "
        " -c 0.9 -gap 100 -g 1 -n 3 -M 0) >> {log}"

rule concat_files:
    output:
        "{dir}/barcodes.clstr"
    input:
        expand("{{dir}}/unique_bc/{sample}.clustered.clstr", sample=indexes)
    shell:
        "cat {input} >> {output}"

rule bowtie2_mapping:
    output:
        bam = "{dir}/mapped.bam"
    input:
        r1_fastq = "{dir}/reads.1.fastq.trimmed.fastq.gz",
        r2_fastq = "{dir}/reads.2.fastq.trimmed.fastq.gz"
    threads: 20
    params:
        reference = paths["bowtie2_reference"]
    log: "{dir}/bowtie2_mapping.log"
    shell:
        " (bowtie2 "
        "    -1 {input.r1_fastq} "
        "    -2 {input.r2_fastq} "
        "    -x {params.reference} "
        "    --maxins 2000 "
        "    -p {threads} | "
        "    samtools view  - "
        "        -@ {threads} "
        "        -bh > {output.bam}) 2> {log}"

rule sort_bam:
    output:
        bam = "{dir}/mapped.sorted.bam"
    input:
        bam = "{dir}/mapped.bam"
    threads: 20
    shell:
        "samtools sort "
        " {input.bam} "
        " -@ {threads} > {output.bam}"

rule tagbam:
    output:
        bam = "{dir}/mapped.sorted.tag.bam"
    input:
        bam = "{dir}/mapped.sorted.bam",
        clstr = "{dir}/barcodes.clstr"
    log: "{dir}/tag_bam.stderr"
    shell:
        "(blr tagbam "
        " {input.bam} "
        " {input.clstr}"
        " {output.bam}"
        " -bc {cluster_tag}) 2> {log} "

rule duplicates_removal:
    output:
        bam = "{dir}/mapped.sorted.tag.rmdup.bam"
    input:
        bam = "{dir}/mapped.sorted.tag.bam"
    log:
        metrics = "{dir}/picard_rmdup_metrics.log",
        stderr = "{dir}/4_rmdup.log"
    shell:
        "(picard MarkDuplicates "
        " I={input.bam} "
        " O={output.bam} "
        " M={log.metrics} "
        " ASSUME_SORT_ORDER=coordinate "
        " REMOVE_DUPLICATES=true "
        " BARCODE_TAG={cluster_tag}) 2> {log.stderr} "

rule duplicates_marking:
    output:
        bam = "{dir}/mapped.sorted.tag.rmdup.mkdup.bam"
    input:
        bam = "{dir}/mapped.sorted.tag.rmdup.bam"
    log:
        metrics = "{dir}/picard_mkdup_metrics.log",
        stderr = "{dir}/4_rmdup.log"
    shell:
        "(picard MarkDuplicates "
        " I={input.bam} "
        " O={output.bam} "
        " M={log.metrics} "
        " ASSUME_SORT_ORDER=coordinate) 2> {log.stderr} "

rule clusterrmdup_and_index:
    output:
        bam = "{dir}/mapped.sorted.tag.rmdup.x2.bam",
        bai = "{dir}/mapped.sorted.tag.rmdup.x2.bam.bai"
    input:
        bam = "{dir}/mapped.sorted.tag.rmdup.mkdup.bam"
    log: "{dir}/4_rmdup.log"
    shell:
        "blr clusterrmdup "
        " {input.bam}"
        " - "
        " -bc {cluster_tag} 2>> {log} | tee {output.bam} | samtools index - {output.bai} "

rule filterclusters:
    output:
        bam = "{dir}/mapped.sorted.tag.rmdup.x2.filt.bam",
        stat1 = "{dir}/cluster_stats/x2.stats.molecules_per_bc",
        stat2 = "{dir}/cluster_stats/x2.stats.molecule_stats"
    input:
        bam = "{dir}/mapped.sorted.tag.rmdup.x2.bam"
    log: "{dir}/4_rmdup.log"
    params:
        stats = "{dir}/cluster_stats/x2.stats"
    shell:
        "(blr filterclusters "
        " -M 260"
        " -s {params.stats} "
        " -bc {cluster_tag} "
        " {input.bam}"
        " {output.bam}) 2>> {log}"
