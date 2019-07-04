"""
Takes a fastq file barcode sequences in the header and writes a barcode fasta file with only unique entries.
"""

import pysam
import logging
import time
from tqdm import tqdm

logger = logging.getLogger(__name__)


def main(args):

    # Generate dict with bc => bc_cluster consensus sequence
    logger.info(f'Starting analysis')

    with open(args.input_clstr, "r") as clstr_file:
        cluster_dict = process_clusters(clstr_file, skip_nonclust=args.skip_nonclust)

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
                read.set_tag('BC', str(bc_id), value_type='Z')
                read.query_name = f'{read.query_name}_BC:Z:{bc_id}'

            out.write(read)

    logger.info(f'Finished')


def process_clusters(input_file, skip_nonclust=False):
    """
    Builds bc => bc_cluster dict (bc_cluster is given as the consensus sequence).
    """

    # For first loop
    seqs_in_cluster = 2 if skip_nonclust else int()

    # Reads cluster file and saves as dict
    cluster_dict = dict()
    cluster_id = int()
    for line in tqdm(input_file, desc="Reading .clstr"):

        # Reports cluster to master dict and start new cluster instance
        if line.startswith('>Cluster'):

            # If non-clustered sequences are to be omitted, removes if only one sequence makes out the cluster
            if skip_nonclust and seqs_in_cluster < 2:
                del cluster_dict[current_key]
            seqs_in_cluster = 0

            cluster_id += 1
            current_value = cluster_id
        else:
            current_key = line.split()[2].lstrip('>').rstrip('...').split(':')[2]
            cluster_dict[current_key] = current_value
            seqs_in_cluster += 1

    return cluster_dict


def add_arguments(parser):
    parser.add_argument("input_mapped_bam", help=".bam file with mapped reads which is to be tagged with barcode id:s.")
    parser.add_argument("input_clstr", help=".clstr file from cdhit clustering.")
    parser.add_argument("output_tagged_bam", help=".bam file with barcode cluster id in the bc tag.")
    parser.add_argument("-s", "--skip_nonclust", action="store_true", help="Does not give cluster ID:s to clusters "
                                                                           "made out by only one sequence.")
