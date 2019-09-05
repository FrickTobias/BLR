"""
Removes barcode tags present at more than -M loci (corresponding to removing barcode tags from reads origin to droplets
which had more than -M molecules in one and the same droplet).
"""

import pysam
import logging

from tqdm import tqdm

logger = logging.getLogger(__name__)


def main(args):
    summary = Summary(args.tags_to_be_removed)
    logger.info(f"Starting")

    # Writes filtered out
    with pysam.AlignmentFile(args.input, "rb") as openin, \
            pysam.AlignmentFile(args.output, "wb", template=openin) as openout:
        logger.info("Writing filtered bam file")
        for read in tqdm(openin.fetch(until_eof=True)):
            summary.tot_reads += 1
            no_mols = fetch_tag(pysam_read=read, no_mols_tag=args.molecule_number_tag, default=0)

            # If barcode is not in all_molecules the barcode does not have enough proximal reads to make a single
            # molecule. If the barcode has more than <max_molecules> molecules, remove it from the read.
            if no_mols > args.max_molecules:
                read, summary = strip_barcode(pysam_read=read, tags_to_be_removed=args.tags_to_be_removed,
                                              summary=summary)

            openout.write(read)

    summary.print_stats()

    logger.info(f"Finished")


def fetch_tag(pysam_read, no_mols_tag, default=None):
    """
    Fetches barcode from a bam file tag, returns None if reads isn't tagged.
    """

    try:
        no_mols = pysam_read.get_tag(no_mols_tag)
    except KeyError:
        no_mols = default

    return no_mols


def strip_barcode(pysam_read, tags_to_be_removed, summary):
    """
    Strips an alignment from its barcode sequence. Keeps information in header but adds FILTERED prior to bc info.
    """

    # Modify header
    header, bc_info = pysam_read.query_name.split("_", maxsplit=1)
    pysam_read.query_name = header + "_" + "FILTERED-" + bc_info

    # Remove tags
    for bam_tag in tags_to_be_removed:
        # Stats
        removed_tag = pysam_read.get_tag(bam_tag)
        summary.removal_dict[bam_tag].add(removed_tag)
        summary.reads_with_removed_tags += 1

        # Strip read from tag
        pysam_read.set_tag(bam_tag, "FILTERED", value_type="Z")

    return pysam_read, summary


class Summary:
    """
    Gathers all stats generated during analysis
    """

    def __init__(self, tags_to_be_removed):

        self.tot_reads = int()
        self.reads_with_removed_tags = int()
        self.removal_dict = dict()
        for bam_tag in tags_to_be_removed:
            self.removal_dict[bam_tag] = set()

    def print_stats(self):
        """
        Prints stats to terminal
        """

        # Read stats
        logger.info(f"Total Reads in file:\t{self.tot_reads:,}")
        for bam_tag, removed_set in self.removal_dict.items():
            logger.info(f"Unique {bam_tag} tags removed: {len(removed_set)}")
        logger.info(f"Reads with barcodes removed:\t{self.reads_with_removed_tags} "
                    f"({self.reads_with_removed_tags / self.tot_reads:.2%})")


def add_arguments(parser):
    parser.add_argument("input", help=".bam file tagged with BX:Z:<int> tags. Needs to be indexed, sorted & have "
                                      "duplicates removed.")
    parser.add_argument("output", help="Output filtered file.")

    parser.add_argument("-M", "--max_molecules", metavar="<INTEGER>", type=int, default=500,
                        help="Maximum number of molecules allowed to keep barcode. DEFAULT: 500")
    parser.add_argument("-mn", "--molecule_number_tag", metavar="<STRING>", type=str, default="MN",
                        help=".bam file tag where number of molecules for that barcode is stored. DEFAULT: MN")
    parser.add_argument("-t", "--tags_to_be_removed", metavar="<STRING>", type=str, nargs="*",
                        default=["BX", "MI", "MN"],
                        help=".bam file tags which are to be removed for filtered reads. DEFAULT: BX MI MN")
