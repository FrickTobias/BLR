"""
Merges a BAM file and cdhit .clster file into an output BAM file which
contains cluster id and barcode sequence under specified BAM tags.
"""

import pysam
import logging
import re
from tqdm import tqdm

logger = logging.getLogger(__name__)


def main(args):

    # Generate dict with bc => bc_cluster consensus sequence
    logger.info(f'Starting analysis')

    with open(args.input_clstr, "r") as clstr_file:
        cluster_dict = process_clusters(clstr_file)

    # Read bam files and translate bc seq to BC cluster ID + write to out
    reads_with_non_clustered_bc = int()
    with pysam.AlignmentFile(args.input_mapped_bam, 'rb') as infile, \
            pysam.AlignmentFile(args.output_tagged_bam, 'wb', template=infile) as out:

        for read in tqdm(infile.fetch(until_eof=True), desc="Reading .bam"):
            read_bc = read.query_name.split()[0].split('_')[-1]

            # Fetch barcode cluster ID based on barcode sequence
            if read_bc not in cluster_dict:
                reads_with_non_clustered_bc += 1
            else:
                bc_id = cluster_dict[read_bc]

                # Stores as string, makes duplicate removal possible. Can do it as integer as well.
                read.set_tag(args.barcode_cluster_tag, str(bc_id), value_type='Z')
                read.query_name = f'{read.query_name}_{args.barcode_cluster_tag}:Z:{bc_id}'
            out.write(read)

    if reads_with_non_clustered_bc:
        logger.info(f'Number of reads not clustered: {reads_with_non_clustered_bc:,}')
    logger.info(f'Finished')


def process_clusters(clstr_file):
    """
    Builds and returns a dictionary of barcode sequences which within the same cluster point to a
    common cluster id number.
    :param clstr_file: open input .clstr file from cdhit
    """

    # Reads cluster file and saves as dict
    cluster_dict = dict()
    cluster_id = int()

    # Assumes lines in format: '0       20nt, >4:1:AACAGTTCTAAATGTGTACA... *'
    # Extracts barcode with letters A,T,C,G of length 18-22 bp between ':' and '...'.
    pattern = re.compile(r":([ATCG]{18,22})...")

    for line in tqdm(clstr_file, desc="Reading .clstr"):

        # Reports cluster to master dict and start new cluster instance
        if line.startswith('>Cluster'):
            cluster_id += 1
        else:
            cluster_sequence = pattern.search(line).group(1)
            cluster_dict[cluster_sequence] = cluster_id

    return cluster_dict


def add_arguments(parser):
    parser.add_argument("input_mapped_bam", metavar="<INPUT_BAM>",
                        help=".bam file with mapped reads which is to be tagged with barcode id:s.")
    parser.add_argument("input_clstr",  metavar="<INPUT_CLSTR>",
                        help=".clstr file from cdhit clustering.")
    parser.add_argument("output_tagged_bam",  metavar="<OUTPUT_BAM>",
                        help=".bam file with barcode cluster id in the bc tag.")
    parser.add_argument("-bc", "--barcode-cluster-tag", metavar="<STRING>", type=str, default="BX",
                        help="Bam file tag where barcode cluster id is stored. 10x genomics longranger output "
                             "uses 'BX' for their error corrected barcodes. DEFAULT: BX")
