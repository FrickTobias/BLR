"""
Transfers SAM tags from query_names to SAM tags. Tags in query_names must follow SAM tag format, e.g. BC:Z:<string>.
"""

import pysam
import logging
import re
from tqdm import tqdm
from collections import Counter

from blr.utils import print_stats

logger = logging.getLogger(__name__)

ALLOWED_BAM_TAG_TYPES = "ABfHiZ"  # From SAM format specs https://samtools.github.io/hts-specs/SAMtags.pdf


def main(args):
    logger.info("Starting analysis")
    summary = Counter()

    # Read SAM/BAM files and transfer barcode information from alignment name to SAM tag
    with pysam.AlignmentFile(args.input, "rb") as infile, \
            pysam.AlignmentFile(args.output, "wb", template=infile) as out:
        for read in tqdm(infile.fetch(until_eof=True), desc="Reading input"):
            for tag in args.tags:

                # Search for tag and return
                match = find_tag(header=read.query_name, bam_tag=tag)

                # If tag string is found, remove from header and set SAM tag value
                if match:
                    full_tag_string = match.group(0)
                    divider = read.query_name[match.start() - 1]
                    read.query_name = read.query_name.replace(divider + full_tag_string, "")
                    read.set_tag(match.group("tag"), match.group("value"), value_type=match.group("type"))
                    summary[f"reads with tag {tag}"] += 1

            out.write(read)

    print_stats(summary, name="stats")
    logger.info("Finished")


def find_tag(header, bam_tag, allowed_value_chars="ATGCN"):
    """
    Finds BAM tags in header and returns regex match objects.
    :param header: strm pysam header
    :param bam_tag: str, SAM tag to search for
    :param allowed_value_chars: str, characters allowed in SAM value
    :return:
    """

    # Build regex pattern strings
    pattern_tag_value = f"[{allowed_value_chars}]+"
    pattern_tag_types = f"[{ALLOWED_BAM_TAG_TYPES}]"

    # Add strings to name match object variables to match.group(<name>)
    regex_bam_tag = add_regex_name(pattern=bam_tag, name="tag")
    regex_allowed_types = add_regex_name(pattern=pattern_tag_types, name="type")
    regex_tag_value = add_regex_name(pattern=pattern_tag_value, name="value")

    # regex pattern search for tag:type:value
    regex_string = r":".join([regex_bam_tag, regex_allowed_types, regex_tag_value])
    return re.search(regex_string, header)


def add_regex_name(pattern, name):
    """
    Formats a string to fit the regex pattern for getting named match objects
    :param pattern: str, regex pattern
    :param name: str
    :return: named regex string
    """
    prefix = "(?P<"
    infix = ">"
    suffix = ")"
    return prefix + name + infix + pattern + suffix


def add_arguments(parser):
    parser.add_argument("input",
                        help="SAM/BAM file with mapped reads which is to be tagged with barcode information. To read "
                             "from stdin use '-'.")

    parser.add_argument("-o", "--output", default="-",
                        help="Write output BAM to file rather then stdout.")
    parser.add_argument("-t", "--tags", default=["BX"], nargs="*", help="List of SAM tags. Default: %(default)s")
