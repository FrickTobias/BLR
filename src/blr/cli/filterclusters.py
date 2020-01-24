"""
Removes barcode tags present at more than -M loci (corresponding to removing barcode tags from reads origin to droplets
which had more than -M molecules in one and the same droplet).
"""

import logging
from collections import Counter
from tqdm import tqdm

from blr.utils import get_bamtag, PySAMIO, print_stats

logger = logging.getLogger(__name__)


def main(args):
    tags_to_remove = [args.barcode_tag, args.molecule_tag, args.number_tag]
    removed_tags = {tag: set() for tag in tags_to_remove}
    summary = Counter()
    logger.info("Starting")

    # Writes filtered out
    with PySAMIO(args.input, args.output, __name__) as (openin, openout):
        for read in tqdm(openin.fetch(until_eof=True), desc="Filtering input", unit="reads"):
            summary["Total reads"] += 1
            no_mols = get_bamtag(pysam_read=read, tag=args.number_tag)

            # If barcode is not in all_molecules the barcode does not have enough proximal reads to make a single
            # molecule. If the barcode has more than <max_molecules> molecules, remove it from the read.
            if no_mols and no_mols > args.max_molecules:
                # Stats
                summary["Removed tags"] += len(tags_to_remove)
                summary["Reads with removed tags"] += 1

                strip_barcode(pysam_read=read, tags_to_be_removed=tags_to_remove, removed_tags=removed_tags)

            openout.write(read)

    summary.update({f"Unique {tag} tags removed": len(removed_tags[tag]) for tag in tags_to_remove})

    logger.info("Finished")

    print_stats(summary, name=__name__)


def strip_barcode(pysam_read, tags_to_be_removed, removed_tags):
    """
    Strips an alignment from its barcode sequence. Keeps information in header but adds FILTERED prior to bc info.
    """

    # Modify header
    pysam_read.query_name = f"{pysam_read.query_name}_FILTERED"

    # Remove tags
    for bam_tag in tags_to_be_removed:
        removed_tags[bam_tag].add(pysam_read.get_tag(bam_tag))
        # Strip read from tag
        pysam_read.set_tag(bam_tag, None, value_type="Z")


def add_arguments(parser):
    parser.add_argument("input",
                        help="SAM/BAM file tagged with barcodes information under the tag specified at "
                             "-b/--barcode-tag. The file needs to be indexed, sorted & have duplicates removed. "
                             "To read from stdin use '-'.")

    parser.add_argument("-o", "--output", default="-",
                        help="Write output BAM to file rather then stdout.")
    parser.add_argument("-b", "--barcode-tag", default="BX",
                        help="SAM tag for storing the error corrected barcode. Default: %(default)s")
    parser.add_argument("-M", "--max_molecules", type=int, default=500,
                        help="Maximum number of molecules allowed to keep barcode. Default: %(default)s")
    parser.add_argument("-m", "--molecule-tag", default="MI",
                        help="SAM tag for storing molecule index specifying a identified molecule for each barcode. "
                             "Default: %(default)s")
    parser.add_argument("-n", "--number_tag", default="MN",
                        help="SAM tag for storing molecule count for a particular barcode. Default: %(default)s")
