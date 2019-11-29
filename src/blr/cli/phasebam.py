"""
Phase BAM files according to phased SNVs
"""

import pysam
import logging
from tqdm import tqdm
from collections import Counter

from blr.utils import print_stats

logger = logging.getLogger(__name__)


def main(args):
    logger.info("Starting analysis")
    summary = Counter()

    # PARSE AND SAVE VCF PHASE INFO

    # Read SAM/BAM files and transfer barcode information from alignment name to SAM tag
    with pysam.AlignmentFile(args.input, "rb") as infile, \
            pysam.AlignmentFile(args.output, "wb", template=infile) as out:
        for read in tqdm(infile.fetch(until_eof=True), desc="Reading input"):

            # ADD PHASING TAGS

            out.write(read)

    print_stats(summary, name="stats")
    logger.info("Finished")


def add_arguments(parser):
    parser.add_argument("input-bam",
                        help="BAM file. To read from stdin use '-'.")
    parser.add_argument("input-vcf",
                        help="Phased VCF file.")

    parser.add_argument("-o", "--output", default="-",
                        help="Write output phased BAM to file rather then stdout.")

