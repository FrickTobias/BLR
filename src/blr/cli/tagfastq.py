"""
Tag FASTQ/FASTA headers with corrected barcodes from separate raw barcode FASTQ with matching headers and starcode
 CLSTR output file with corrected barcodes.
"""

import logging
import sys
import dnaio
from tqdm import tqdm
from collections import namedtuple

logger = logging.getLogger(__name__)


def main(args):
    logger.info(f'Starting')

    # Get the corrected barcodes and create a dictionary pointing each raw barcode to its canonical sequence.

    with open(args.barcodes_corrected, "r") as reader:
        corrected_barcodes = parse_corrected_barcodes(reader)

    # Get the raw barcodes and create a dictionary pointing each read header to the raw barcode sequence.
    with dnaio.open(args.barcodes_raw, mode="r") as reader:
        header_barcodes = parser_raw_barcodes(reader)

    in_interleaved = True if not args.input2 else False
    logger.info(f"Input detected as {'interleaved fastq.' if in_interleaved else 'paired fastq.'}")

    # If no output1 is given output is sent to stdout
    if not args.output1:
        logger.info(f"Writing output to stdout.")
        args.output1 = sys.stdout.buffer
        args.output2 = None

    out_interleaved = True if not args.output2 else False
    logger.info(f"Output detected as {'interleaved fastq.' if out_interleaved else 'paired fastq.'}")

    reads_skipped = list()
    # Parse input FASTA/FASTQ for read1 and read2 and write output
    with dnaio.open(args.input1, file2=args.input2, interleaved=in_interleaved, mode="r") as reader, \
            dnaio.open(args.output1, file2=args.output2, interleaved=out_interleaved, mode="w") as writer:
        for read1, read2 in tqdm(reader, desc='Read pairs processed'):
            # Header parsing
            name_and_pos_r1, read_and_index_r1 = read1.name.split(maxsplit=1)
            name_and_pos_r2, read_and_index_r2 = read2.name.split(maxsplit=1)
    
            assert name_and_pos_r1 == name_and_pos_r2
    
            try:
                raw_barcode_seq = header_barcodes[name_and_pos_r1]
            except KeyError:
                reads_skipped.append(name_and_pos_r1)
                continue
            
            corr_barcode = corrected_barcodes[raw_barcode_seq]

            # Make new string with barcode sequence and id to add to headers.
            new_text = f"{corr_barcode.seq}_{args.barcode_cluster_tag}:Z:{corr_barcode.id}"

            # Save header to read instances
            read1.name = f"{name_and_pos_r1}_{new_text} {read_and_index_r1}"
            read2.name = f"{name_and_pos_r2}_{new_text} {read_and_index_r2}"
    
            # Write to out
            writer.write(read1, read2)

    for header in reads_skipped:
        logger.info(f"Barcode for read-pair {header} not found. Skipped read-pair.")

    logger.info(f'Finished')


def parse_corrected_barcodes(open_file):
    """
    Parse starcode cluster output and return a dictionary with raw sequences pointing to a
    corrected canonical sequence
    :param open_file: starcode tabular output file.
    :return: dict: raw sequences pointing to a nametuple containe corrected canonical sequence and id.
    """
    target = namedtuple("Cluster", ['seq', 'id'])
    corrected_barcodes = dict()
    for id, cluster in tqdm(enumerate(open_file.readlines(), start=1), desc="Corrected barcodes processed"):
        canonical_seq, _, cluster_seqs = cluster.strip().split("\t", maxsplit=3)
        cluster_target = target(seq=canonical_seq, id=id)
        corrected_barcodes.update({raw_seq: cluster_target for raw_seq in cluster_seqs.split(",")})
    return corrected_barcodes


def parser_raw_barcodes(open_file):
    """
    Parse FASTA/FASTQ containing barcode sequences and return a dictionary with entry headers pointing to a
    raw barcodes sequence.
    :param open_file: dnaio odject.
    :return: dict: entry headers pointing to a raw barcodes sequence
    """
    header_barcodes = dict()
    for barcode in tqdm(open_file, desc="Raw barcodes processed"):
        header, _ = barcode.name.split(maxsplit=1)
        header_barcodes[header] = barcode.sequence
    return header_barcodes


def add_arguments(parser):
    parser.add_argument(
        "input1", metavar='<INPUT1>',
        help="Input FASTQ/FASTA file. Assumes to contain read1 if given with second input file. If only "
             "input1 is given, input is assumed to be an interleaved. If reading from stdin"
             "is requested use '-' as a placeholder.")
    parser.add_argument(
        "input2", nargs='?', metavar='<INPUT2>',
        help="Input  FASTQ/FASTA for read2 for paired-end read. Leave empty if using interleaved.")
    parser.add_argument(
        "-o1", "--output1", default=None, metavar='<OUTPUT1>',
        help="Output FASTQ/FASTA file name for read1. If not specified the result is written to stdout"
             " as interleaved. If output1 given but not output2, output will be written as "
             "interleaved to output1.")
    parser.add_argument(
        "-o2", "--output2", default=None, metavar='<OUTPUT2>',
        help="Output FASTQ/FASTA name for read2. If not specified but -o1 given the result is written"
             " as interleaved to o1.")
    parser.add_argument(
        "-rb", "--barcodes-raw", required=True, metavar='<RAW BARCODES>',
        help="FASTQ/FASTA for raw barcodes.")
    parser.add_argument(
        "-cb", "--barcodes-corrected", required=True, metavar='<CORRECTED BARCODES>',
        help="FASTQ/FASTA for error corrected barcodes. Currently accepts output from starcode clustering with "
             "'--print-clusters' enabled.")
    parser.add_argument(
        "-bc", "--barcode-cluster-tag", metavar="<STRING>", type=str, default="BX",
        help="BAM file tag where barcode cluster id is stored. 10x genomics longranger output "
             "uses 'BX' for their error corrected barcodes. DEFAULT: BX")
