
variants = "variants.reference.vcf" if config["reference_variants"] else "variants.called.vcf"
bamfile_basename = "mapped.sorted.tag.mkdup.bcmerge.mol.filt.BQSR" if config["BQSR"] else "mapped.sorted.tag.mkdup.bcmerge.mol.filt"

rule hapcut2_extracthairs:
    output:
        unlinked = bamfile_basename + ".unlinked.txt"
    input:
        bam = bamfile_basename + ".bam",
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
        linked = bamfile_basename + ".linked.txt"
    input:
        bam = bamfile_basename + ".bam",
        bai = bamfile_basename + ".bam.bai",
        vcf = variants,
        unlinked = bamfile_basename + ".unlinked.txt"
    log: "hapcut2_linkfragments.log"
    shell:
         "LinkFragments.py"
         " --bam {input.bam}"
         " -v {input.vcf}"
         " --fragments {input.unlinked}"
         " --out {output.linked} &> {log}"


rule hapcut2_phasing:
    output:
        phase = bamfile_basename + ".phase",
        phased_vcf = bamfile_basename + ".phase.phased.VCF"
    input:
        linked = bamfile_basename + ".linked.txt",
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
         vcf1 = bamfile_basename + ".phase.phased.VCF",
         vcf2 = "ground_truth.phased.vcf"
    shell:
         "blr calculate_haplotype_statistics"
         " -v1 {input.vcf1}"
         " -v2 {input.vcf2}"
         " > {output.stats}"
