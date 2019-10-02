"""
Transfers cluster id and barcode sequence from BAM file header to BAM tags in new output BAM.
"""

import pysam
import logging
from tqdm import tqdm

logger = logging.getLogger(__name__)


def main(args):

    # Generate dict with bc => bc_cluster consensus sequence
    logger.info(f'Starting analysis')

    # Read bam files and translate bc seq to BC cluster ID + write to out
    with pysam.AlignmentFile(args.input_mapped_bam, 'rb') as infile, \
            pysam.AlignmentFile(args.output_tagged_bam, 'wb', template=infile) as out:

        for read in tqdm(infile.fetch(until_eof=True), desc="Reading .bam"):
            read_header, cluster_seq, cluster_tag = read.query_name.split()[0].split('_')
            cluster_id = cluster_tag.split(":")[-1]

            read.set_tag(args.barcode_cluster_tag, cluster_id, value_type='Z')
            out.write(read)

    logger.info(f'Finished')


def add_arguments(parser):
    parser.add_argument("input_mapped_bam", metavar="<INPUT_BAM>",
                        help=".bam file with mapped reads which is to be tagged with barcode id:s.")
    parser.add_argument("output_tagged_bam",  metavar="<OUTPUT_BAM>",
                        help=".bam file with barcode cluster id in the bc tag.")
    parser.add_argument("-bc", "--barcode-cluster-tag", metavar="<STRING>", type=str, default="BX",
                        help="Bam file tag where barcode cluster id is stored. 10x genomics longranger output "
                             "uses 'BX' for their error corrected barcodes. DEFAULT: BX")
