import pandas as pd
variants = "variants.reference.vcf" if config["reference_variants"] else "variants.called.vcf"

def chr_from_file(file):
    with open(file) as f:
        return [line.strip() for line in f]

# Get list of chromosomes for phasing analysis.
chromosomes = [f"chr{nr}" for nr in range(22)] if not config["chromosomes"] else chr_from_file(config["chromosomes"])


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


rule bzgip_and_index_vcf:
    output:
        vcf = "mapped.sorted.tag.mkdup.bcmerge.mol.filt.phase.phased.vcf.gz",
        index = "mapped.sorted.tag.mkdup.bcmerge.mol.filt.phase.phased.vcf.gz.tbi"
    input:
        vcf = "mapped.sorted.tag.mkdup.bcmerge.mol.filt.phase.phased.VCF"
    shell:
        "bgzip -c {input.vcf} > {output.vcf} &&"
        " tabix -p vcf {output.vcf}"


rule split_vcf:
    output: expand("chromosome_phased_vcf/{chromosome}.phased.vcf", chromosome=chromosomes)
    input: "mapped.sorted.tag.mkdup.bcmerge.mol.filt.phase.phased.vcf.gz"
    run:
        for chromosome in chromosomes:
            shell(
                "tabix {input} {chromosome} > {output}"
            )

rule hapcut2_stats:
    output:
        stats = "chromosome_phased_vcf/{chromosome}.phasing_stats.txt"
    input:
         vcf1 = "chromosome_phased_vcf/{chromosome}.phased.vcf"
    params:
          vcf2 = config["phasing_ground_truth"]
    shell:
         "calculate_haplotype_statistics.py"
         " -v1 {input.vcf1}"
         " -v2 {params.vcf2}"
         " > {output.stats}"

rule aggregate_stats:
    output: "phasing_stats.csv"
    input: expand("chromosome_phased_vcf/{chromosome}.phasing_stats.txt", chromosome=chromosomes)
    run:
        results = list()
        for file in input:
            row = {"chromosome": file.split("/")[1].split(".")[0]}
            with open(file, "r") as f:
                for l in f:
                    if ":" not in l:
                        continue
                    row[l.split(":")[0]] = l.split(":")[1].strip()
            results.append(row)
        df = pd.DataFrame(results).set_index("chromosome")
        df.to_csv(output[0])
