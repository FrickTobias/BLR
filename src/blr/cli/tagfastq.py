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
from itertools import islice
from tqdm import tqdm

logger = logging.getLogger(__name__)


def main(args):
    logger.info("Starting")

    # Get the corrected barcodes and create a dictionary pointing each raw barcode to its
    # canonical sequence.
    with open(args.corrected_barcodes, "r") as reader:
        corrected_barcodes = parse_corrected_barcodes(reader)

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
    separator = args.sep
    # Parse input FASTA/FASTQ for read1 and read2, uncorrected barcodes and write output
    with dnaio.open(args.input1, file2=args.input2, interleaved=in_interleaved, mode="r",
                    fileformat="fastq") as reader, \
            dnaio.open(args.output1, file2=args.output2, interleaved=out_interleaved, mode="w",
                       fileformat="fastq") as writer, \
            BarcodeReader(args.uncorrected_barcodes) as uncorrected_barcode_reader:

        for read1, read2 in tqdm(reader, desc="Read pairs processed"):
            # Header parsing
            name_and_pos_r1, read_and_index_r1 = read1.name.split(maxsplit=1)
            name_and_pos_r2, read_and_index_r2 = read2.name.split(maxsplit=1)

            uncorrected_barcode_seq = uncorrected_barcode_reader.get_barcode(name_and_pos_r1)

            # Check if barcode was found and update header with barcode info.
            if uncorrected_barcode_seq:
                corrected_barcode_seq = corrected_barcodes[uncorrected_barcode_seq]

                raw_barcode_id = f"{args.sequence_tag}:Z:{uncorrected_barcode_seq}"
                corr_barcode_id = f"{args.barcode_tag}:Z:{corrected_barcode_seq}"

                # Create new name with barcode information.
                new_name = separator.join([name_and_pos_r1, raw_barcode_id, corr_barcode_id])

                # Save header to read instances
                read1.name = " ".join([new_name, read_and_index_r1])
                read2.name = " ".join([new_name, read_and_index_r2])
            else:
                reads_missing_barcode += 1

            # Write to out
            writer.write(read1, read2)

    logger.info(f"Read-pairs missing barcodes: {reads_missing_barcode}")

    logger.info("Finished")


def parse_corrected_barcodes(open_file):
    """
    Parse starcode cluster output and return a dictionary with raw sequences pointing to a
    corrected canonical sequence
    :param open_file: starcode tabular output file.
    :return: dict: raw sequences pointing to a corrected canonical sequence.
    """
    corrected_barcodes = dict()
    for cluster in tqdm(open_file, desc="Clusters processed"):
        canonical_seq, _, cluster_seqs = cluster.strip().split("\t", maxsplit=3)
        corrected_barcodes.update({raw_seq: canonical_seq for raw_seq in cluster_seqs.split(",")})
    return corrected_barcodes


class BarcodeReader:
    def __init__(self, filename):
        self._cache = dict()
        self._file = dnaio.open(filename, mode="r")
        self.barcodes = self.parse()

    def parse(self):
        for barcode in tqdm(self._file, desc="Uncorrected barcodes processed"):
            read_name, _ = barcode.name.split(maxsplit=1)
            yield (read_name, barcode.sequence)

    def get_barcode(self, read_name, maxiter=10):
        if read_name in self._cache:
            return self._cache.pop(read_name)

        for barcode_read_name, barcode_sequence in islice(self.barcodes, maxiter):
            # If read_name in next pair then parser lines are synced --> drop cache.
            if read_name == barcode_read_name:
                self._cache.clear()
                return barcode_sequence

            self._cache[barcode_read_name] = barcode_sequence
        return None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self._file.close()


def add_arguments(parser):
    parser.add_argument(
        "uncorrected_barcodes",
        help="FASTQ/FASTA for uncorrected barcodes.")
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
        "-b", "--barcode-tag", default="BX",
        help="SAM tag for storing the error corrected barcode. Default: %(default)s")
    parser.add_argument(
        "-s", "--sequence-tag", default="RX",
        help="SAM tag for storing the uncorrected barcode sequence. Default: %(default)s")
    parser.add_argument(
        "--sep", default="_",
        help="Character used as separator for storing SAM tags in the FASTQ/FASTA header. Default: %(default)s"
    )
