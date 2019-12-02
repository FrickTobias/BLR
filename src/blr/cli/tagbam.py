"""
Transfers SAM tags from query headers to SAM tags. Currently tags in header must follow SAM tag format, e.g.
BC:Z:<SEQUENCE>.
"""

import pysam
import logging
import re
from tqdm import tqdm
from collections import Counter

from blr.utils import print_stats

logger = logging.getLogger(__name__)

ALLOWED_SAM_TAG_TYPES = "ABfHiZ"  # From SAM format specs https://samtools.github.io/hts-specs/SAMtags.pdf


def main(args):
    # Can't be at top since function are defined later
    FUNCTION_DICT = {
        "rx-bx": rx_bx,
        "ema": ema
    }

    logger.info("Starting analysis")
    summary = Counter()
    adjust_header_and_set_tag = FUNCTION_DICT[args.format]

    # Read SAM/BAM files and transfer barcode information from alignment name to SAM tag
    with pysam.AlignmentFile(args.input, "rb") as infile, \
            pysam.AlignmentFile(args.output, "wb", template=infile) as out:
        for read in tqdm(infile.fetch(until_eof=True), desc="Reading input", unit=" reads"):

            adjust_header_and_set_tag(read, args.separator)

            for tag in tags:
                read.query_name.split(args.separator, maxsplit=1)[0]
                if not args.only_remove:
                    read.set_tag(tag, match.group("value"), value_type=match.group("type"))
                summary[f"Reads with tag {tag}"] += 1
            summary["Total read count"] += 1

            out.write(read)

    print_stats(summary, name=__name__)
    logger.info("Finished")


def rx_bx(read, separator, first_tag, second_tag):
    """
    Adjusts headers for any format
    :param read:
    :param separator:
    :return:
    """

    try:
        name, rx, bx = read.query_name.split(separator)
    except ValueError:
        return

    read.query_name = name
    for tag in [rx,bx]
        tag, tag_type, val = tag.split(":")
        read.set_tag(tag, val, value_type=tag_type)


def ema(read)

    read.query_name = "".join(query_name.rsplit(":", 1))


def add_arguments(parser):
    parser.add_argument("input",
                        help="BAM file with SAM tag info in header. To read from stdin use '-'.")

    parser.add_argument("-o", "--output", default="-",
                        help="Write output BAM to file rather then stdout.")
    parser.add_argument("-f", "--format", default="two_tags", choices=["two_tags", "ema"],
                        help="Specify what tag search function to use for finding tags. Default: %(default)s")
    parser.add_argument("-s", "--separator")
    parser.add_argument("--tag-list", )
