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
        fastq="{dir}/reads.fastq.trimmed.tag.sorted.itlvd.fastq"
    input:
        r1="{dir}/reads.1.fastq.trimmed.tag.fastq.gz",
        r2="{dir}/reads.2.fastq.trimmed.tag.fastq.gz"
    shell:
         """
         blr sortfastq {input.r1} {input.r2} {output.fastq} --interleaved_output
         """

rule megahit:
    output:
        fa="{dir}/contigs.fa"
    input:
        fq="{dir}/reads.fastq.trimmed.tag.sorted.itlvd.fastq"
    threads: 4
    params:
        tmpdir="$TMPDIR/megahit"
    shell:
        """
        rm -rf {params.tmpdir}
        megahit --12 {input.fq} -o {params.tmpdir} -t {threads}
        mv {params.tmpdir}/final.contigs.fa {output.fa}
        rm -r {params.tmpdir}
        """

rule bwa_index:
    output:
        index=expand("{{dir}}/contigs.fa.{suffix}", suffix=["amb","ann","bwt","pac","sa"])
    input:
        fa="{dir}/contigs.fa"
    shell:
        """
        bwa index {input.fa}
        """

rule bwa_mem:
    output:
        bam="{dir}/reads.read_cloud_preprocessings.bam"
    input:
        index=expand("{{dir}}/contigs.fa.{suffix}", suffix=["amb","ann","bwt","pac","sa"]),
        fastq="{dir}/reads.fastq.trimmed.tag.sorted.itlvd.fastq"
    params:
        index_base="{dir}/contigs.fa"
    shell:
        """
        bwa mem -C -p {input.fastq} {params.index_base} | samtools sort -o {output.bam} -
        """