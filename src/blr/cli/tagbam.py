"""
Transfers cluster id and barcode sequence from BAM file header to BAM tags in new output BAM.
"""

import pysam
import logging
from tqdm import tqdm

logger = logging.getLogger(__name__)


def main(args):

    # Generate dict with bc => bc_cluster consensus sequence
    logger.info("Starting analysis")

    alignments_missing_bc = 0

    # Read bam files and translate bc seq to BC cluster ID + write to out
    with pysam.AlignmentFile(args.input_mapped_bam, 'rb') as infile, \
            pysam.AlignmentFile(args.output_tagged_bam, 'wb', template=infile) as out:

        for read in tqdm(infile.fetch(until_eof=True), desc="Reading BAM"):
            try:
                cluster_tag = read.query_name.split()[0].split('_')[-1]
            except ValueError:
                alignments_missing_bc += 1
                cluster_tag = None

            if cluster_tag:
                cluster_id = cluster_tag.split(":")[-1]

                read.set_tag(args.barcode_tag, cluster_id, value_type='Z')

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
