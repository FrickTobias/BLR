HAPCUT2=config["hapcut2"]

rule HAPCUT2_extractHAIRS:
    output:
        unlinked = "{dir}/mapped.sorted.tag.rmdup.x2.filt.unlinked"
    input:
        bam = "{dir}/mapped.sorted.tag.rmdup.x2.filt.bam",
        vcf = "{dir}/reference.vcf"
    log: "{dir}/hapcut2_extracthairs.log"
    shell:
         "{HAPCUT2}/build/extractHAIRS"
         " --10X 1"
         " --bam {input.bam}"
         " --VCF {input.vcf}"
         " --out {output.unlinked} 2> {log}"

rule HAPCUT2_LinkFragments:
    output:
        linked = "{dir}/mapped.sorted.tag.rmdup.x2.filt.linked"
    input:
        bam = "{dir}/mapped.sorted.tag.rmdup.x2.filt.bam",
        vcf = "{dir}/reference.vcf",
        unlinked = "{dir}/mapped.sorted.tag.rmdup.x2.filt.unlinked"
    log: "{dir}/hapcut2_linkfragments.log"
    shell:
         "python {HAPCUT2}/utilities/LinkFragments.py"
         " --bam {input.bam}"
         " -v {input.vcf}"
         " --fragments {input.unlinked}"
         " --out {output.linked} &> {log}"

rule HAPCUT2_phasing:
    output:
        phase = "{dir}/mapped.sorted.tag.rmdup.x2.filt.phase",
        phased_vcf = "{dir}/mapped.sorted.tag.rmdup.x2.filt.phase.phased.vcf"
    input:
        linked = "{dir}/mapped.sorted.tag.rmdup.x2.filt.linked",
        vcf = "{dir}/reference.vcf"
    log: "{dir}/hapcut2_phasing.log"
    shell:
         "{HAPCUT2}/build/HAPCUT2"
         " --nf 1"
         " --fragments {input.linked}"
         " --vcf {input.vcf}"
         " --out {output.phase}"
         " --outvcf 1 2> {log}"
