rule tagfastq:
    """Tags each read in a fastq file with the consensus barcode sequence"""
    output:
        r1="{dir}/reads.1.fastq.trimmed.tag.fastq.gz",
        r2="{dir}/reads.2.fastq.trimmed.tag.fastq.gz"
    input:
        r1="{dir}/reads.1.fastq.trimmed.fastq.gz",
        r2="{dir}/reads.2.fastq.trimmed.fastq.gz",
        clstr="{dir}/BC.NNN.clstr"
    params:
        r1="{dir}/reads.1.fastq.trimmed.tag.fastq",
        r2="{dir}/reads.2.fastq.trimmed.tag.fastq"
    shell:
        """
        python python\ scripts/tag_fastq.py {input.r1} {input.r2} {input.clstr} {params.r1} {params.r2}
        gzip {params.r1}
        gzip {params.r2}
        """