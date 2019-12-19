
variants = "variants.reference.vcf" if config["reference_variants"] else "variants.called.vcf"

rule hapcut2_extracthairs:
    output:
        unlinked = "mapped.sorted.tag.mkdup.bcmerge.mol.filt.unlinked.txt"
    input:
        bam = "mapped.sorted.tag.mkdup.bcmerge.mol.filt.bam",
        vcf = variants
    log: "hapcut2_extracthairs.log"
    shell:
         "extractHAIRS"
         " --10X 1"
         " --bam {input.bam}"
         " --VCF {input.vcf}"
         " --out {output.unlinked} 2> {log}"


rule hapcut2_linkfragments:
    output:
        linked = "mapped.sorted.tag.mkdup.bcmerge.mol.filt.linked.txt"
    input:
        bam = "mapped.sorted.tag.mkdup.bcmerge.mol.filt.bam",
        vcf = variants,
        unlinked = "mapped.sorted.tag.mkdup.bcmerge.mol.filt.unlinked.txt"
    log: "hapcut2_linkfragments.log"
    shell:
         "LinkFragments.py"
         " --bam {input.bam}"
         " -v {input.vcf}"
         " --fragments {input.unlinked}"
         " --out {output.linked} &> {log}"


rule hapcut2_phasing:
    output:
        phase =      "mapped.sorted.tag.mkdup.bcmerge.mol.filt.phase",
        phased_vcf = "mapped.sorted.tag.mkdup.bcmerge.mol.filt.phase.phased.VCF"
    input:
        linked = "mapped.sorted.tag.mkdup.bcmerge.mol.filt.linked.txt",
        vcf = variants
    log: "hapcut2_phasing.log"
    shell:
         "hapcut2"
         " --nf 1"
         " --fragments {input.linked}"
         " --vcf {input.vcf}"
         " --out {output.phase}"
         " --outvcf 1 2> {log}"

rule symlink_reference_phased:
    output: "ground_truth.phased.vcf"
    shell: "ln -s {config[phasing_ground_truth]} {output}"

rule hapcut2_stats:
    output:
        stats = "phasing_stats.txt"
    input:
         vcf1 = "mapped.sorted.tag.mkdup.bcmerge.mol.filt.phase.phased.VCF",
         vcf2 = "ground_truth.phased.vcf"
    shell:
         "blr calculate_haplotype_statistics"
         " -v1 {input.vcf1}"
         " -v2 {input.vcf2}"
         " > {output.stats}"
