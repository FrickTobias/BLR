"""
Adds SAM tags to all reads in a BAM/SAM/CRAM file.

For RG tags it also modifies the header appropriately, assuming:
  - File only contains one library
  - File contains Illumina short reads

Usage example:
  blr addbamtag --tag RG --value 4 --type Z non-tagged.bam tagged.bam

"""

import pysam
import logging
from collections import Counter
from tqdm import tqdm

from blr import utils

logger = logging.getLogger(__name__)

RG_TEMPLATE = {
    'ID': '',
    'LB': 'lib1',
    'PU': 'unit1',
    'SM': '20',
    'PL': 'ILLUMINA'}


def main(args):
    logger.info(f"Starting analysis")
    tag_value = args.value if not args.value == "None" else None

    # Go through all reads in input fq
    with pysam.AlignmentFile(args.input, "rb") as reader:
        header = add_tag_to_header(args.input, args.tag, args.value) if args.tag == "RG" else reader.header
        with pysam.AlignmentFile(args.output, "wb", header=header) as output:
            for read in tqdm(reader):
                read.set_tag(args.tag, tag_value, value_type=args.type)
                output.write(read)

    logger.info(f"Finished")


def add_tag_to_header(bamfile, tag, value):
    if tag != "RG":
        logger.critical(f"Trying to create RG header for non-RG tag addition (tag: {tag}, value: {value}).")
        exit("Exiting")
    with pysam.AlignmentFile(bamfile, "rb") as openin:
        rg_header_dict = RG_TEMPLATE.copy()
        rg_header_dict["ID"] = value

        rg_header = list()
        rg_header.append(rg_header_dict)

        new_header = openin.header.to_dict()
        new_header["RG"] = rg_header

    return pysam.AlignmentHeader.from_dict(new_header)


def add_arguments(parser):
    parser.add_argument("input",
                        help="SAM/BAM file.")
    parser.add_argument("-o", "--output", default="-",
                        help="Write output BAM to file rather then stdout.")
    parser.add_argument("--tag", required=True, help="SAM tag which will be written. Required.")
    parser.add_argument("--type", default="Z", help="SAM tag type character. %(default)s.")
    parser.add_argument("--value", required=True,
                        help="Value which will be added to SAM tag. If set to 'None' removes tag entirely instead. "
                             "Required.")
