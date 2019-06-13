rule tagfastq:
    """Tags each read with a cluster number from barcode clustering"""
    output:
        r1="{dir}/reads.1.fastq.trimmed.tag.fastq.gz",
        r2="{dir}/reads.2.fastq.trimmed.tag.fastq.gz"
    input:
        r1="{dir}/reads.1.fastq.trimmed.fastq.gz",
        r2="{dir}/reads.2.fastq.trimmed.fastq.gz",
        clstr="{dir}/BC.NNN.clstr"
    shell:
        """
        blr tagfastq {input.r1} {input.r2} {input.clstr} {output.r1} {output.r2}
        """

rule sortfastq:
    output:
        r1="{dir}/reads.1.fastq.trimmed.tag.sorted.fastq.gz",
        r2="{dir}/reads.2.fastq.trimmed.tag.sorted.fastq.gz"
    input:
        r1="{dir}/reads.1.fastq.trimmed.tag.fastq.gz",
        r2="{dir}/reads.2.fastq.trimmed.tag.fastq.gz"
    shell:
         """
         blr sortfastq {input.r1} {input.r2} {output.r1} --output2 {output.r2}
         """

