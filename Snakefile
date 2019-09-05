from snakemake.utils import validate
import itertools
import os

configfile: "config.yaml"
validate(config, "config.schema.yaml")

# Create list of files to be created in cdhitprep
indexes = sorted(["".join(tup) for tup in itertools.product("ATCG", repeat=config["index_nucleotides"])]) \
                if config["index_nucleotides"] > 0 else ["all"]

# Import rules for trimming fastq files.
include: "rules/trim.smk"

rule compress:
    output: "{dir}/{sample}.fastq.gz"
    input: "{dir}/{sample}.fastq"
    shell:
        "pigz < {input} > {output}"

# If the number of index nucleotide is 0 only on file will be created.
if config["index_nucleotides"] == 0:
    rule cdhitprep_no_index:
        # Create fasta containing aggregates barcode sequences from fastq file headers.
        output:
            "{dir}/unique_bc/all.fa"
        input:
            r1_fastq = "{dir}/trimmed-c.1.fastq.gz"
        log:
            stdout = "{dir}/cdhit_prep.stdout",
            stderr = "{dir}/cdhit_prep.stderr"
        shell:
            "blr cdhitprep"
            " {input.r1_fastq}"
            " {output}"
            " -f 0 > {log.stdout} 2> {log.stderr}"
else:
    rule cdhitprep:
        # Create fasta containing aggregates barcode sequences from fastq file headers.
        output:
            expand("{{dir}}/unique_bc/{sample}.fa", sample=indexes)
        input:
            r1_fastq = "{dir}/trimmed-c.1.fastq.gz"
        params:
            dir = "{dir}/unique_bc/"
        log:
            stdout = "{dir}/cdhit_prep.stdout",
            stderr = "{dir}/cdhit_prep.stderr"
        shell:
            "blr cdhitprep"
            " {input.r1_fastq}"
            " {params.dir}"
            " -i {config[index_nucleotides]}"
            " -f 0 > {log.stdout} 2> {log.stderr}"


rule barcode_clustering:
    # Barcode clustering using cd-hit-454
    input:
       "{dir}/unique_bc/{sample}.fa"
    output:
        "{dir}/unique_bc/{sample}.clustered",
        "{dir}/unique_bc/{sample}.clustered.clstr"
    threads: 20
    log: "{dir}/unique_bc/{sample}.clustering.log"
    params:
        prefix= lambda wildcards,output: os.path.splitext(output[1])[0]
    shell:
        "cd-hit-454"
        " -i {input}"
        " -o {params.prefix}"
        " -T {threads}"
        " -c 0.9"
        " -gap 100"
        " -g 1"
        " -n 3"
        " -M 0 >> {log}"

rule concat_files:
    # Concatenate all the .clstr files into one single file.
    output:
        "{dir}/barcodes.clstr"
    input:
        expand("{{dir}}/unique_bc/{sample}.clustered.clstr", sample=indexes)
    shell:
        "cat {input} > {output}"

rule bowtie2_mapping:
    # Mapping of trimmed fastq to reference using bowtie2 and sorting output using samtools.
    output:
        bam = "{dir}/mapped.sorted.bam"
    input:
        r1_fastq = "{dir}/trimmed-c.1.fastq.gz",
        r2_fastq = "{dir}/trimmed-c.2.fastq.gz"
    threads: 20
    params:
        reference = config["bowtie2_reference"]
    log: "{dir}/bowtie2_mapping.log"
    shell:
        "bowtie2"
        " -1 {input.r1_fastq}"
        " -2 {input.r2_fastq}"
        " -x {params.reference}"
        " --maxins 2000"
        " -p {threads} 2> {log} |"
        "samtools sort -"
        " -@ {threads}"
        " -O BAM > {output.bam}"

rule tagbam:
    # Add barcode information to bam file using custom script
    output:
        bam = "{dir}/mapped.sorted.tag.bam"
    input:
        bam = "{dir}/mapped.sorted.bam",
        clstr = "{dir}/barcodes.clstr"
    log: "{dir}/tag_bam.stderr"
    shell:
        "blr tagbam"
        " {input.bam}"
        " {input.clstr}"
        " {output.bam}"
        " -bc {config[cluster_tag]} 2> {log}"

rule buildmolecules:
    # Groups reads into molecules depending on their genomic position and barcode
    output:
        bam = "{dir}/mapped.sorted.tag.mol.bam"
    input:
        bam = "{dir}/mapped.sorted.tag.bam"
    log: "{dir}/buildmolecules.log"
    shell:
        "blr buildmolecules"
        " {input.bam}"
        " {output.bam}"
        " -bc {config[cluster_tag]}"
        " 2> {log}"

rule duplicates_marking:
    # Mark duplicates between barcode clusters using picard
    output:
        bam = "{dir}/mapped.sorted.tag.mol.mkdup.bam"
    input:
        bam = "{dir}/mapped.sorted.tag.mol.bam"
    log:
        metrics = "{dir}/picard_mkdup_metrics.log",
        stderr = "{dir}/4_rmdup.log"
    params:
        picard_command = config["picard_command"],
        heap_space=config["heap_space"]
    shell:
        "{params.picard_command} -Xms{params.heap_space}g MarkDuplicates"
        " I={input.bam}"
        " O={output.bam}"
        " M={log.metrics}"
        " ASSUME_SORT_ORDER=coordinate 2> {log.stderr}"

rule clusterrmdup_and_index:
    # Removes cluster duplicates and indexes output
    output:
        bam = "{dir}/mapped.sorted.tag.mol.mkdup.bcmerge.bam",
        bai = "{dir}/mapped.sorted.tag.mol.mkdup.bcmerge.bam.bai",
        merges = "{dir}/barcode-merges.csv"
    input:
        bam = "{dir}/mapped.sorted.tag.mol.mkdup.bam"
    log: "{dir}/4_rmdup.log"
    shell:
        "blr clusterrmdup"
        " {input.bam}"
        " -"
        " {output.merges}"
        " -bc {config[cluster_tag]} 2>> {log} |"
        "tee {output.bam} |"
        "samtools index - {output.bai}"

rule filterclusters:
    # Filter clusters based on parameters
    output:
        bam = "{dir}/mapped.sorted.tag.mol.mkdup.bcmerge.filt.bam",
    input:
        bam = "{dir}/mapped.sorted.tag.mol.mkdup.bcmerge.bam"
    log: "{dir}/filterclusters.log"
    shell:
        "blr filterclusters"
        " {input.bam}"
        " {output.bam}"
        " -M 260"
        " -t {config[cluster_tag]} 2> {log}"

rule bam_to_fastq:
    # Convert final bam file to fastq files for read 1 and 2
    output:
        r1_fastq = "{dir}/reads.1.final.fastq",
        r2_fastq = "{dir}/reads.2.final.fastq"
    input:
        bam = "{dir}/mapped.sorted.tag.mol.mkdup.bcmerge.filt.bam"
    log: "{dir}/picard_samtofastq.log"
    params:
        picard_command = config["picard_command"],
        heap_space=config["heap_space"]
    shell:
        "{params.picard_command} -Xms{params.heap_space}g SamToFastq"
        " I={input.bam}"
        " FASTQ={output.r1_fastq}"
        " SECOND_END_FASTQ={output.r2_fastq} 2>> {log}"
