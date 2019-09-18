"""
Extract barcode sequences by moving 20 bp from 5' end of the read to the header,
separated by an underline. Only accepts paired-end input (interleaved or not).
Script can handle gzipped files (.gz).
"""
import logging
import sys
import dnaio
from tqdm import tqdm

logger = logging.getLogger(__name__)


def main(args):
    logger.info(f'Starting')

    input_interleaved = True if not args.input2 else False
    logger.info(f"Input detected as {'interleaved fastq.' if input_interleaved else 'paired fastq.'}")

    # If no output1 is given output is sent to stdout
    if not args.output1:
        logger.info(f"Writing output to stdout.")
        args.output1 = sys.stdout.buffer
        args.output2 = None

    output_interleaved = True if not args.output2 else False
    logger.info(f"Output detected as {'interleaved fastq.' if output_interleaved else 'paired fastq.'}")

    reader = dnaio.open(args.input1, file2=args.input2, interleaved=input_interleaved, mode="r", fileformat="fastq")
    writer = dnaio.open(args.output1, file2=args.output2, interleaved=output_interleaved, mode="w", fileformat="fastq")
    for read1, read2 in tqdm(reader, desc='Read pairs processed'):
        # Adjusting for BC
        bc_seq = read1.sequence[:20]
        read1.sequence = read1.sequence[20:]
        read1.qualities = read1.qualities[20:]

        # Header parsing
        name_and_pos_r1, read_and_index_r1 = read1.name.split(maxsplit=1)
        name_and_pos_r2, read_and_index_r2 = read2.name.split(maxsplit=1)

        # Save header to read instances
        read1.name = name_and_pos_r1 + '_' + bc_seq + ' ' + read_and_index_r1
        read2.name = name_and_pos_r2 + '_' + bc_seq + ' ' + read_and_index_r2

        # Write to out
        writer.write(read1, read2)

    reader.close()
    writer.close()
    logger.info(f'Finished')


def add_arguments(parser):
    parser.add_argument(
        "input1", metavar='<INPUT FASTQ1>',
        help="Input .fastq file. Assumes to contain read1 if given with second input file. If only "
             "input1 is given, input is assumed to be an interleaved fastq file. If reading from stdin"
             "is requested use '-' as a placeholder.")
    parser.add_argument(
        "input2", nargs='?', metavar='<INPUT FASTQ2>',
        help="Input .fastq file for read2 for paired-end read. Leave empty if using interleaved fastq.")
    parser.add_argument(
        "-o1", "--output1", default=None, metavar='<OUTPUT FASTQ1>',
        help="Output .fastq file name for read1. If not specified the result is written to stdout"
             " as interleaved fastq. If output1 given but not output2, output will be written as "
             "interleaved fastq to output1.")
    parser.add_argument(
        "-o2", "--output2", default=None, metavar='<OUTPUT FASTQ2>',
        help="Output .fastq file name for read2. If not specified but -o1 given the result is written"
             " as interleaved fastq to o1.")
