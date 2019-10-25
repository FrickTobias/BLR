"""
Transfers barcode sequence informatino from BAM file header to BAM tags in new output BAM.
"""

import pysam
import logging
from tqdm import tqdm

logger = logging.getLogger(__name__)


def main(args):

    # Generate dict with bc => bc_cluster consensus sequence
    logger.info("Starting analysis")

    alignments_missing_bc = 0

    # Read BAM files and transfer barcode information from alignment name to SAM tag
    with pysam.AlignmentFile(args.input_mapped_bam, 'rb') as infile, \
            pysam.AlignmentFile(args.output_tagged_bam, 'wb', template=infile) as out:

        for read in tqdm(infile.fetch(until_eof=True), desc="Reading BAM"):
            try:
                name, raw_barcode_tag, corr_barcode_tag = read.query_name.split()[0].split('_')
            except ValueError:
                alignments_missing_bc += 1
                name, raw_barcode_tag, corr_barcode_tag = read.query_name.split()[0], None, None

            # Rename aligment to exclude barcode information.
            read.query_name = name

            if corr_barcode_tag:
                assert corr_barcode_tag.startswith(f"{args.barcode_tag}:Z:")
                corr_barcode = corr_barcode_tag.split(":")[-1]
                read.set_tag(args.barcode_tag, corr_barcode, value_type='Z')

            if raw_barcode_tag:
                assert raw_barcode_tag.startswith(f"{args.sequence_tag}:Z:")
                raw_barcode = raw_barcode_tag.split(":")[-1]
                read.set_tag(args.sequence_tag, raw_barcode, value_type='Z')

            out.write(read)

    logger.info(f"Alignments missing barcodes: {alignments_missing_bc}")
    logger.info("Finished")


def add_arguments(parser):
    parser.add_argument("input_mapped_bam",
                        help="BAM file with mapped reads which is to be tagged with barcode ids.")
    parser.add_argument("output_tagged_bam",
                        help="BAM file with barcode cluster id in the bc tag.")
    parser.add_argument("-b", "--barcode-tag", default="BX",
                        help="SAM tag for storing the error corrected barcode. Default: %(default)s")
    parser.add_argument("-s", "--sequence-tag", default="RX",
                        help="SAM tag for storing the raw barcode sequence. Default: %(default)s")
