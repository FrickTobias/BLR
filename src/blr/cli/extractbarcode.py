"""
Extract barcode sequences by moving 20 bp from 5' end of the read to the header,
separated by an underline.
"""
import logging

import blr.utils as BLR

logger = logging.getLogger(__name__)


def main(args):
    logger.info('Starting')
    progress = BLR.ProgressReporter('Read pairs processed', 1000000)
    generator = BLR.FileReader(args.r1, args.r2)
    with open(args.out_r1, 'w') as openr1, open(args.out_r2, 'w') as openr2:
        for read1, read2 in generator.fastqPairedReader():

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
            openr1.write(read1.fastq_string())
            openr2.write(read2.fastq_string())

            # Progress reporting
            progress.update()

    generator.close()
    logger.info('Finished')


def add_arguments(parser):
    parser.add_argument("r1", help="Read 1 fastq file")
    parser.add_argument("r2", help="Read 2 fastq file")
    parser.add_argument("out_r1", help="Read 1 output")
    parser.add_argument("out_r2", help="Read 2 output")
