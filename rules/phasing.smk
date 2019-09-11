HAPCUT2=config["hapcut2"]

rule hapcut2_extracthairs:
    output:
        unlinked = "{dir}/mapped.sorted.tag.mkdup.bcmerge.mol.filt.unlinked.txt"
    input:
        bam = "{dir}/mapped.sorted.tag.mkdup.bcmerge.mol.filt.bam",
        vcf = "{dir}/reference.vcf"
    log: "{dir}/hapcut2_extracthairs.log"
    shell:
         "{HAPCUT2}/build/extractHAIRS"
         " --10X 1"
         " --bam {input.bam}"
         " --VCF {input.vcf}"
         " --out {output.unlinked} 2> {log}"

rule hapcut2_linkfragments:
    output:
        linked = "{dir}/mapped.sorted.tag.mkdup.bcmerge.mol.filt.linked.txt"
    input:
        bam = "{dir}/mapped.sorted.tag.mkdup.bcmerge.mol.filt.bam",
        vcf = "{dir}/reference.vcf",
        unlinked = "{dir}/mapped.sorted.tag.mkdup.bcmerge.mol.filt.unlinked.txt"
    log: "{dir}/hapcut2_linkfragments.log"
    shell:
         "python {HAPCUT2}/utilities/LinkFragments.py"
         " --bam {input.bam}"
         " -v {input.vcf}"
         " --fragments {input.unlinked}"
         " --out {output.linked} &> {log}"

rule hapcut2_phasing:
    output:
        phase =      "{dir}/mapped.sorted.tag.mkdup.bcmerge.mol.filt.phase",
        phased_vcf = "{dir}/mapped.sorted.tag.mkdup.bcmerge.mol.filt.phase.phased.VCF"
    input:
        linked = "{dir}/mapped.sorted.tag.mkdup.bcmerge.mol.filt.linked.txt",
        vcf = "{dir}/reference.vcf"
    log: "{dir}/hapcut2_phasing.log"
    shell:
         "{HAPCUT2}/build/HAPCUT2"
         " --nf 1"
         " --fragments {input.linked}"
         " --vcf {input.vcf}"
         " --out {output.phase}"
         " --outvcf 1 2> {log}"

rule hapcut2_stats:
    output:
        stats = "{dir}/phasing_stats.txt"
    input:
         vcf1 = "{dir}/mapped.sorted.tag.mkdup.bcmerge.mol.filt.phase.phased.VCF",
    params:
          vcf2 = config["phasing_ground_truth"]
    shell:
         "python {HAPCUT2}/utilities/calculate_haplotype_statistics.py"
         " -v1 {input.vcf1}"
         " -v2 {params.vcf2}"
         " > {output.stats}"