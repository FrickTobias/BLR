"""
Phase BAM files according to phased SNVs
"""

import pysam
import logging
import vcfpy
from tqdm import tqdm
from collections import Counter, OrderedDict

from blr.utils import print_stats, get_bamtag

logger = logging.getLogger(__name__)


def main(args):
    logger.info("Starting analysis")

    # Save hetSNV info
    phased_snv_dict = get_phased_snvs(args.input_vcf)

    # Phase reads & corresponding molecules at hetSNV sites
    phasing_summary = Counter()
    molecule_phasing_dict, phased_single_reads = phase_molecules(args.input_bam, args.molecule_tag, phased_snv_dict,
                                                                 phasing_summary)
    print_stats(phasing_summary, name=f"{__name__} - Run details")

    # Write output (setting phasing)
    output_summary = Counter()
    with pysam.AlignmentFile(args.input_bam, "rb") as infile, \
            pysam.AlignmentFile(args.output, "wb", template=infile) as out:
        for read in tqdm(infile.fetch(until_eof=True), desc="Writing output", unit=" reads"):
            output_summary["Total reads"] += 1

            # Set phasing tag
            molecule = get_bamtag(read, args.molecule_tag)
            if molecule in molecule_phasing_dict:
                haplotype = max(molecule_phasing_dict[molecule])
                read.set_tag(args.haplotype_tag, haplotype)
                output_summary["Phased reads in molecules"] += 1
            elif read.query_name in phased_single_reads:
                haplotype = phased_single_reads[read.query_name]
                read.set_tag(args.haplotype_tag, haplotype)
                output_summary["Phased single reads"] += 1

            out.write(read)

    print_stats(output_summary, name=__name__)
    logger.info("Finished")


def get_phased_snvs(vcf_file):
    """
    Reads VCF files and saves all phased SNVs.
    :param vcf_file: file, vcf format
    :return: dict, dict[chrom][pos][nucleotide] = phase_info
    """
    vcf_reader = vcfpy.Reader.from_path(vcf_file)
    phased_snv_dict = dict()
    for record in tqdm(vcf_reader, desc="Reading VCF file", unit=" records processed"):
        chrom = record.CHROM
        ref_pos = record.POS

        # Save phased SNVs
        for call in record.calls:
            if not call.is_phased:
                continue

            # Setup dict for position
            if chrom not in phased_snv_dict:
                phased_snv_dict[chrom] = OrderedDict()

            # Get call (1 or 0) and translate to nucleotide
            h1 = call.gt_bases[call.gt_alleles[0]]
            h2 = call.gt_bases[call.gt_alleles[1]]

            # Save phasing info
            if ref_pos not in phased_snv_dict[chrom]:
                phased_snv_dict[chrom][ref_pos] = dict()
            phased_snv_dict[chrom][ref_pos][h1] = "h1"  # TODO: Change this to HapCUT2 flags for phase info
            phased_snv_dict[chrom][ref_pos][h2] = "h2"  # TODO: Change this to HapCUT2 flags for phase info

    return phased_snv_dict


def phase_molecules(bam_file, molecule_tag, phased_snv_dict, summary):
    """
    Goes through a SAM file and matches reads at variants to one of the two haplotypes. Then sets the reads molecule
    phase accordingly.
    :param bam_file: file, SAM format
    :param molecule_tag: str, SAM tag
    :param phased_snv_dict: dict, dict[chrom][pos][nucleotide] = phase_info
    :param summary: Counter instance, from collections
    :return: dict, dict[molecule] = phase_info
    """
    molecule_phasing_dict = dict()
    phased_single_reads = dict()
    prev_read_pos = (str(), int(), int())
    with pysam.AlignmentFile(bam_file, "rb") as infile:
        for read in tqdm(infile.fetch(until_eof=True), desc="Phasing reads at hetSNV positions", unit=" reads"):

            # Init criterias for read
            if skip_read(read, phased_snv_dict, summary):
                continue

            # Update variant positions if new read start/end
            if prev_read_pos != (read.reference_name, read.reference_start, read.reference_end):
                variants_at_read, phased_snv_dict = get_phased_variants(read, phased_snv_dict)
                prev_read_pos = (read.reference_name, read.reference_start, read.reference_end)

            # Link molecule to phasing information
            if len(variants_at_read) >= 1:

                # Phase read to the phase with most support. Returns None if read does not match called alleles.
                haplotype = phase_read(read, variants_at_read, summary)
                if not haplotype:
                    continue

                # Increment support for molecule belonging to haplotype
                molecule = get_bamtag(read, molecule_tag)
                if molecule:
                    if molecule not in molecule_phasing_dict:
                        molecule_phasing_dict[molecule] = Counter()
                    molecule_phasing_dict[molecule][haplotype] += 1
                    summary["Phased reads with molecule info"] += 1
                # If no molecule for read, only phase read pair
                else:
                    phased_single_reads[read.query_name] = haplotype
                    summary["Phased reads without molecule info"] += 1

    return molecule_phasing_dict, phased_single_reads


def skip_read(read, phased_snv_dict, summary):
    if read.reference_start is not None and read.reference_end is None:
        summary["Read with ref start but no ref end"] += 1
        return True
    if read.reference_name not in phased_snv_dict:
        summary["Reads without phased SNV in chr"] += 1
        return True
    # if first variant is after read, continue to next read
    if read.reference_end < first_item(phased_snv_dict[read.reference_name]):
        return True


def first_item(ordered):
    return next(iter(ordered))


def get_phased_variants(read, phased_snv_dict):
    """
    Finds all called hetSNVs withing read.reference_start and read.reference_end. Removes variants downstream of
    read.reference_start from phased_snv_dict.

    :param read: pysam read alignment
    :param phased_snv_dict: dict, dict[chrom][pos][nucleotide] = phase_info, must be sorted within chromosomes
    :return: dict[pos][nucleotide] = phase_info, updated phased_snv_dict
    """
    variants_at_read = dict()
    chrom = read.reference_name

    for var_pos, alleles in phased_snv_dict[chrom].copy().items():

        # if first variant is before ref start, remove var from dict
        if read.reference_start > var_pos:
            del phased_snv_dict[chrom][var_pos]
            continue
        elif read.reference_start <= var_pos and read.reference_end >= var_pos:
            variants_at_read[var_pos] = alleles
    return variants_at_read, phased_snv_dict


def phase_read(read, variants_at_read, summary):
    """
    Assigns a read a haplotype from the variants at that read
    :param read: pysam read alignment
    :param variants_at_read: dict, dict[ref_pos][phase_info] = allele_seq
    :param summary: Counter instance, from collections
    :return: str, phase_info for the variant with most support.
    """

    # Translate to read.seq positions
    pos_translations = translate_positions(read=read, ref_pos_list=list(variants_at_read.keys()))
    read_haplotype = Counter()

    for ref_pos, seq_pos in pos_translations.items():
        # Using read.query_alignment_sequence compensates for insertions (rather than read.seq)
        nt = read.query_alignment_sequence[seq_pos]

        if nt in variants_at_read[ref_pos]:
            haplotype = variants_at_read[ref_pos][nt]
            read_haplotype[haplotype] += 1

    # Read did not fit any of the called haplotypes
    if len(read_haplotype) == 0:
        summary["Reads not fitting called alleles"] += 1
        return

    return max(read_haplotype)


def translate_positions(read, ref_pos_list):
    """
    Translates a reference position to the nucleotide position in an read.seq string, where deletions are accounted
    for.

    Discrepancy problem between Ref pos/read.seq pos. ref start =/= read start & deletions.

    Ref        =======================
    Ref pos     1 3     9
    Read        |-- . . ---> (. = Deletion)
    read pos    0 2     3

    :param ref_pos_list: list with ints, position in reference to be translated
    :param read: pysam read alignment
    :return: dict, dict[ref_pos] = seq_pos
    """

    ref_pos_iterator = iter(ref_pos_list)
    ref_pos = next(ref_pos_iterator)
    pos_translations = dict()

    # Loop adjusting for deletion (extra bases in ref / missing bases in read)
    seq_pos = int()
    for block in read.get_blocks():
        block_start, block_stop = block
        while True:

            # Pos upstream of block, add block len to position
            if block_stop < ref_pos:
                seq_pos += block_stop - block_start
                break
            # Pos between blocks
            elif block_start > ref_pos:
                pos_translations[ref_pos] = None  # TODO: Change this to whatever is used to describe deletions

            # Pos found
            else:
                block_pos = ref_pos - block_start
                pos_translations[ref_pos] = seq_pos + block_pos - 1

            # When last variant, return (even if not last alignment block)
            try:
                ref_pos = next(ref_pos_iterator)
            except StopIteration:
                return pos_translations


def add_arguments(parser):
    parser.add_argument("input_bam",
                        help="BAM file. To read from stdin use '-'. Must be sorted.")
    parser.add_argument("input_vcf",
                        help="Phased VCF file. Must be sorted.")

    parser.add_argument("-o", "--output", default="-",
                        help="Write output phased BAM to file rather then stdout.")
    parser.add_argument("--molecule-tag", default="MI", help="Molecule SAM tag. Default: %(default)s.")
    parser.add_argument("--haplotype-tag", default="HP", help="Haplotype SAM tag. Default: %(default)s")
