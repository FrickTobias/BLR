"""
Tag FASTQ/FASTA headers with corrected barcodes from separate raw barcode FASTQ with matching
headers and starcode CLSTR output file with corrected barcodes.

ABOUT:

First the raw barcodes FASTQ is parser to get a dictionary:

    raw_barcodes[<HEADER>] = <RAW_BARCODE>

Then the corrected barcodes CLSTR file from starcode are parsed to get a dictionary:

    corrected_barcodes[<RAW_BARCODE>] = <CORRECTED_BARCODE>

The for each read-pair in the input FASTQ(s) the corrected barcode is recovered and used to tag
the read by including it in the header.

    <HEADER> ==> <RAW_BARCODE> ==> <CORRECTED_BARCODE>
"""

import logging
import sys
import dnaio
from tqdm import tqdm
from collections import namedtuple

logger = logging.getLogger(__name__)


def main(args):
    logger.info("Starting")

    # Get the corrected barcodes and create a dictionary pointing each raw barcode to its
    # canonical sequence.
    with open(args.corrected_barcodes, "r") as reader:
        corrected_barcodes = parse_corrected_barcodes(reader)

    # Get the raw barcodes and create a dictionary pointing each read header to the raw
    # barcode sequence.
    with dnaio.open(args.raw_barcodes, mode="r") as reader:
        raw_barcodes = parse_raw_barcodes(reader)

    in_interleaved = not args.input2
    logger.info(f"Input detected as {'interleaved' if in_interleaved else 'paired'} FASTQ.")

    # If no output1 is given output is sent to stdout
    if not args.output1:
        logger.info("Writing output to stdout.")
        args.output1 = sys.stdout.buffer
        args.output2 = None

    out_interleaved = not args.output2
    logger.info(f"Output detected as {'interleaved' if out_interleaved else 'paired'} FASTQ.")

    reads_missing_barcode = 0
    # Parse input FASTA/FASTQ for read1 and read2 and write output
    with dnaio.open(args.input1, file2=args.input2, interleaved=in_interleaved, mode="r") as reader, \
            dnaio.open(args.output1, file2=args.output2, interleaved=out_interleaved, mode="w") as writer:
        for read1, read2 in tqdm(reader, desc="Read pairs processed"):
            # Header parsing
            name_and_pos_r1, read_and_index_r1 = read1.name.split(maxsplit=1)
            name_and_pos_r2, read_and_index_r2 = read2.name.split(maxsplit=1)

            try:
                raw_barcode_seq = raw_barcodes[name_and_pos_r1]
            except KeyError:
                raw_barcode_seq = None
                reads_missing_barcode += 1

            if raw_barcode_seq:
                corr_barcode = corrected_barcodes[raw_barcode_seq]

                # Make new string with barcode sequence and id to add to headers.
                new_text = f"{corr_barcode.seq}_{args.barcode_cluster_tag}:Z:{corr_barcode.index}"

                # Save header to read instances
                read1.name = f"{name_and_pos_r1}_{new_text} {read_and_index_r1}"
                read2.name = f"{name_and_pos_r2}_{new_text} {read_and_index_r2}"

            # Write to out
            writer.write(read1, read2)

    logger.info(f"Read-pairs missing barcodes: {reads_missing_barcode}")

    logger.info("Finished")


def parse_corrected_barcodes(open_file):
    """
    Parse starcode cluster output and return a dictionary with raw sequences pointing to a
    corrected canonical sequence
    :param open_file: starcode tabular output file.
    :return: dict: raw sequences pointing to a nametuple containe corrected canonical
    sequence and index.
    """
    target = namedtuple("Cluster", ['seq', 'index'])
    corrected_barcodes = dict()
    for index, cluster in tqdm(enumerate(open_file.readlines(), start=1), desc="Clusters processed"):
        canonical_seq, _, cluster_seqs = cluster.strip().split("\t", maxsplit=3)
        cluster_target = target(seq=canonical_seq, index=index)
        corrected_barcodes.update({raw_seq: cluster_target for raw_seq in cluster_seqs.split(",")})
    return corrected_barcodes


def parse_raw_barcodes(open_file):
    """
    Parse FASTA/FASTQ containing barcode sequences and return a dictionary with entry headers
    pointing to a raw barcodes sequence.
    :param open_file: dnaio odject.
    :return: dict: entry headers pointing to a raw barcodes sequence
    """
    raw_barcodes = dict()
    for barcode in tqdm(open_file, desc="Raw barcodes processed"):
        header, _ = barcode.name.split(maxsplit=1)
        raw_barcodes[header] = barcode.sequence
    return raw_barcodes


def add_arguments(parser):
    parser.add_argument(
        "raw_barcodes",
        help="FASTQ/FASTA for raw barcodes.")
    parser.add_argument(
        "corrected_barcodes",
        help="FASTQ/FASTA for error corrected barcodes. Currently accepts output from starcode "
             "clustering with '--print-clusters' enabled.")
    parser.add_argument(
        "input1",
        help="Input FASTQ/FASTA file. Assumes to contain read1 if given with second input file. "
             "If only input1 is given, input is assumed to be an interleaved. If reading from stdin"
             "is requested use '-' as a placeholder.")
    parser.add_argument(
        "input2", nargs='?',
        help="Input  FASTQ/FASTA for read2 for paired-end read. Leave empty if using interleaved.")
    parser.add_argument(
        "--output1", "--o1",
        help="Output FASTQ/FASTA file name for read1. If not specified the result is written to "
             "stdout as interleaved. If output1 given but not output2, output will be written as "
             "interleaved to output1.")
    parser.add_argument(
        "--output2", "--o2",
        help="Output FASTQ/FASTA name for read2. If not specified but --o1/--output1 given the "
             "result is written as interleaved .")
    parser.add_argument(
        "--barcode-cluster-tag", "--bc", default="BX",
        help="BAM file tag where barcode cluster id is stored. 10x genomics longranger output "
             "uses 'BX' for their error corrected barcodes. Default: %(default)s")
