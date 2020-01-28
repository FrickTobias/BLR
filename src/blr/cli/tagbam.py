"""
Strips headers from tags and depending on mode, set the appropriate SAM tag.
"""

import logging
from tqdm import tqdm
from collections import Counter

from blr.utils import print_stats, PySAMIO

logger = logging.getLogger(__name__)


def main(args):
    # Can't be at top since function are defined later
    function_dict = {
        "bowtie2": mode_samtags_underline_separation,
        "ema": mode_ema,
        "bwa": mode_samtags_underline_separation,
        "minimap2": mode_samtags_underline_separation
    }

    logger.info("Starting analysis")
    summary = Counter()

    # Read SAM/BAM files and transfer barcode information from alignment name to SAM tag
    with PySAMIO(args.input, args.output, __name__) as (infile, outfile):
        mapper = infile.header.to_dict()["PG"][0]["PN"]
        processing_function = function_dict[mapper]

        for read in tqdm(infile.fetch(until_eof=True), desc="Processing reads", unit=" reads"):
            # Strips header from tag and depending on script mode, possibly sets SAM tag
            summary["Total reads"] += 1
            processing_function(read, summary)
            outfile.write(read)

    print_stats(summary, name=__name__)
    logger.info("Finished")


def mode_samtags_underline_separation(read, summary):
    """
    Trims header from tags and sets SAM tags according to values found in header.
    Assumes format: @header_<tag>:<type>:<seq> (can be numerous tags). Constrictions are: Header includes SAM tags
    separated by "_".
    :param read: pysam read alignment
    :param summary: Collections's Counter object
    """

    # Strip header
    header = read.query_name.split("_")
    read.query_name = header[0]

    # Set SAM tags
    for tag in header[1:]:
        tag, tag_type, val = tag.split(":")
        read.set_tag(tag, val, value_type=tag_type)
        summary[f"Reads with tag {tag}"] += 1


def mode_ema(read, *unused):
    """
    Trims header from barcode sequences.
    Assumes format @header:and:more...:header:<seq>. Constrictions: There must be exactly 9 elements separated by ":"
    :param read: pysam read alignment
    """

    # Strip header
    read.query_name = read.query_name.rsplit(":", 1)[0]


def add_arguments(parser):
    parser.add_argument("input",
                        help="BAM file with SAM tag info in header. To read from stdin use '-'.")

    parser.add_argument("-o", "--output", default="-",
                        help="Write output BAM to file rather then stdout.")
