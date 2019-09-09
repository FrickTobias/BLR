HAPCUT2=config["hapcut2"]

rule HAPCUT2_extractHAIRS:
    output:
        unlinked = "{dir}/mapped.sorted.tag.mkdup.bcmerge.filt.unlinked"
    input:
        bam = "{dir}/mapped.sorted.tag.mkdup.bcmerge.filt.bam",
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
        linked = "{dir}/mapped.sorted.tag.mkdup.bcmerge.filt.linked"
    input:
        bam = "{dir}/mapped.sorted.tag.mkdup.bcmerge.filt.bam",
        vcf = "{dir}/reference.vcf",
        unlinked = "{dir}/mapped.sorted.tag.mkdup.bcmerge.filt.unlinked"
    log: "{dir}/hapcut2_linkfragments.log"
    shell:
         "python {HAPCUT2}/utilities/LinkFragments.py"
         " --bam {input.bam}"
         " -v {input.vcf}"
         " --fragments {input.unlinked}"
         " --out {output.linked} &> {log}"

rule HAPCUT2_phasing:
    output:
        phase =      "{dir}/mapped.sorted.tag.mkdup.bcmerge.filt.phase",
        phased_vcf = "{dir}/mapped.sorted.tag.mkdup.bcmerge.filt.phase.phased.VCF"
    input:
        linked = "{dir}/mapped.sorted.tag.mkdup.bcmerge.filt.linked",
        vcf = "{dir}/reference.vcf"
    log: "{dir}/hapcut2_phasing.log"
    params:
        phase = "{dir}/mapped.sorted.tag.mkdup.bcmerge.filt.phase"
    shell:
         "{HAPCUT2}/build/HAPCUT2"
         " --nf 1"
         " --fragments {input.linked}"
         " --vcf {input.vcf}"
         " --out {params.phase}"
         " --outvcf 1 2> {log}"

rule HapCUT2_stats:
    output:
        stats = "{dir}/phasing_stats.txt"
    input:
         vcf1 = "{dir}/mapped.sorted.tag.mkdup.bcmerge.filt.phase.phased.VCF",
    params:
          vcf2 = config["phasing_ground_truth"]
    shell:
         "python {HAPCUT2}/utilities/calculate_haplotype_statistics.py"
         " -v1 {input.vcf1}"
         " -v2 {params.vcf2}"
         " > {output.stats}"