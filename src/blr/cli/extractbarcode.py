"""
Extract barcode sequences by moving 20 bp from 5' end of the read to the header,
separated by an underline.
"""
import logging
import sys

import blr.utils as BLR

logger = logging.getLogger(__name__)


def main(args):
    logger.info(f'Starting')

    progress = BLR.ProgressReporter('Read pairs processed', 1000000)

    # Check if paired-end files given or not.
    input1, input2 = (None, None)
    if len(args.input) == 2:
        input1, input2 = args.input
    elif len(args.input) == 1:
        input1 = args.input[0]
    else:
        sys.exit(f'EXITING: Files given as input ({len(args.input)}) does not match the required 1 or 2.')

    generator = BLR.FileReader(input1, filehandle2=input2, gzipped=args.gzipped, interleaved=args.interleaved)

    read1_name = None
#    for read in generator.fastqReader():
    for read1, read2 in generator.reader():
        # Adjusting for BC
        bc_seq = read1.seq[:20]
        read1.seq = read1.seq[20:]
        read1.qual = read1.qual[20:]

        # Header parsing
        name_and_pos_r1, read_and_index_r1 = read1.header.split(maxsplit=1)
        name_and_pos_r2, read_and_index_r2 = read2.header.split(maxsplit=1)

        # Save header to read instances
        read1.header = name_and_pos_r1 + '_' + bc_seq + ' ' + read_and_index_r1
        read2.header = name_and_pos_r2 + '_' + bc_seq + ' ' + read_and_index_r2

        # Write to out
        sys.stdout.write(read1.fastq_string())
        sys.stdout.write(read2.fastq_string())

        # Progress reporting
        progress.update()

    # with open(args.out_r1, 'w') as openr1, open(args.out_r2, 'w') as openr2:
    #     for read1, read2 in generator.fastqPairedReader():
    #
    #         # Adjusting for BC
    #         bc_seq = read1.seq[:20]
    #         read1.seq = read1.seq[20:]
    #         read1.qual = read1.qual[20:]
    #
    #         # Header parsing
    #         name_and_pos_r1, read_and_index_r1 = read1.header.split(maxsplit=1)
    #         name_and_pos_r2, read_and_index_r2 = read2.header.split(maxsplit=1)
    #
    #         # Save header to read instances
    #         read1.header = name_and_pos_r1 + '_' + bc_seq + ' ' + read_and_index_r1
    #         read2.header = name_and_pos_r2 + '_' + bc_seq + ' ' + read_and_index_r2
    #
    #         # Write to out
    #         openr1.write(read1.fastq_string())
    #         openr2.write(read2.fastq_string())
    #
    #         # Progress reporting
    #         progress.update()

    generator.close()
    logger.info(f'Finished')


def add_arguments(parser):
    # Files
    parser.add_argument("input", nargs='*',
                        help="Input paired fastq for read1 and read2. Use '-' if using stdin.")
    parser.add_argument("-o1", default=None,
                        help="Output file name for read1. If not specified the result is written to stdout as "
                             "interleaved fastq.")
    parser.add_argument("-o2", default=None,
                        help="Output file name for read2. If not specified the result is written to stdout as "
                             "interleaved fastq.")
    # Options
    parser.add_argument("--interleaved", default=False, action='store_true',
                        help="Interleaved fastq")
    parser.add_argument("--gzipped", default=False, action='store_true',
                        help="Input is gzipped. Required if stdin used as input. For other files this may "
                             "be detected from file extension.")
    #parser.add_argument("r2", help="Read 2 fastq file")
    #parser.add_argument("out_r1", help="Read 1 output")
    #parser.add_argument("out_r2", help="Read 2 output")
