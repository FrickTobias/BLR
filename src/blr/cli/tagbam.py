"""
Takes a fastq file barcode sequences in the header and writes a barcode fasta file with only unique entries.
"""

import pysam

import blr.utils as BLR


def main(args):

    # Generate dict with bc => bc_cluster consensus sequence
    BLR.report_progress("Starting analysis")
    clstr_generator = BLR.FileReader(args.input_clstr)
    cluster_dict = ProcessClusters(clstr_generator.fileReader())
    clstr_generator.close()

    # Read bam files and translate bc seq to BC cluster ID + write to out
    progress = BLR.ProgressReporter('Reads processed', 1000000)
    infile = pysam.AlignmentFile(args.input_mapped_bam, 'rb')
    out = pysam.AlignmentFile(args.output_tagged_bam, 'wb', template=infile)
    reads_with_non_clustered_bc = int()
    for read in infile.fetch(until_eof=True):
        read_bc = read.query_name.split()[0].split('_')[-1]

        # Fetch barcode cluster ID based on barcode sequence
        if not read_bc in cluster_dict:
            reads_with_non_clustered_bc += 1
        else:
            bc_id = cluster_dict[read_bc]
            read.set_tag('BC', str(bc_id), value_type='Z')  # Stores as string, makes duplicate removal possible. Can do it as integer as well.
            read.query_name = (read.query_name + '_BC:Z:' + str(bc_id))

        out.write(read)
        progress.update()

    infile.close()
    out.close()
    BLR.report_progress('Finished')

def ProcessClusters(openInfile):
    """
    Builds bc => bc_cluster dict (bc_cluster is given as the consensus sequence).
    """

    # For first loop
    if args.skip_nonclust: seqs_in_cluster = 2

    # Reads cluster file and saves as dict
    cluster_dict = dict()
    cluster_ID = int()
    for line in openInfile:

        # Reports cluster to master dict and start new cluster instance
        if line.startswith('>'):

            # If non-clustered sequences are to be omitted, removes if only one sequence makes out the cluster
            if args.skip_nonclust and seqs_in_cluster < 2:
                del cluster_dict[current_key]
            seqs_in_cluster = 0

            cluster_ID += 1
            current_value = cluster_ID
        else:
            current_key = line.split()[2].lstrip('>').rstrip('...').split(':')[2]
            cluster_dict[current_key] = current_value
            seqs_in_cluster +=1

    return(cluster_dict)

def add_arguments(parser):
    parser.add_argument("input_mapped_bam", help=".bam file with mapped reads which is to be tagged with barcode id:s.")
    parser.add_argument("input_clstr", help=".clstr file from cdhit clustering.")
    parser.add_argument("output_tagged_bam", help=".bam file with barcode cluster id in the bc tag.")
    parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                       "Not recommended due to different function "
                                                                       "names in python 2 and 3.")
    parser.add_argument("-s", "--skip_nonclust", action="store_true", help="Does not give cluster ID:s to clusters "
                                                                           "made out by only one sequence.")