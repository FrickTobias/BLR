"""
Extract barcode sequences by moving 20 bp from 5' end of the read to the header,
separated by an underline.
"""
import logging
import sys
import dnaio

import blr.utils as BLR

logger = logging.getLogger(__name__)


def input_files(inputs, interleaved):
    """
    Checks input arguments to see if they are sound and assigns them to input1 and input2
    :param inputs: list of strings. Input files or paths given.
    :param interleaved: boolean. If input is interleaved or not.
    :return:
    """
    input1, input2 = (None, None)
    if len(inputs) == 0:
        sys.exit(f'EXITING: No inputs given!')
    elif len(inputs) > 2:
        sys.exit(f'EXITING: Too many files given as input {" ".join(inputs)}')
    elif len(inputs) == 2:
        input1, input2 = inputs
        if interleaved:
            interleaved = False
        logger.info(f'Extracting barocodes from fastq {input1} and {input2}')
    elif len(inputs) == 1 and interleaved:
        input1 = inputs[0]
        logger.info(f'Extracting barocodes from interleaved fastq {input1}')
    else:
        sys.exit(f'EXITING: Files given as input ({inputs}) does not match the required 1 or 2.')

    return input1, input2, interleaved


def output_files(output1, output2):
    """
    Checks output arguments given and if file should be written as interleaved or not.
    :param output1: string. Output1 argument
    :param output2: string. Output2 argument
    :return: output1, output2, interleaved.
    """
    interleaved = bool()
    if output1 and output2:
        interleaved = False
        logger.info(f'Writing paired fastq to files {output1} and {output2}')
    elif not output1 and not output2:
        output1 = sys.stdout.buffer
        interleaved = True
        logger.info(f'Writing interleaved fastq to stdout.')
    elif output1 and not output2:
        interleaved = True
        logger.info(f'Writing interleaved fastq to file {output1}.')

    return output1, output2, interleaved


def main(args):
    logger.info(f'Starting')

    progress = BLR.ProgressReporter('Read pairs processed', 1000000)

    inputfile1, inputfile2, input_interleaved = input_files(args.input, args.interleaved)

    reader = dnaio.open(inputfile1, file2=inputfile2, interleaved=input_interleaved, mode="r", fileformat="fastq")

    outputfile1, outputfile2, output_interleaved = output_files(args.output1, args.output2)

    writer = dnaio.open(outputfile1, file2=outputfile2, interleaved=output_interleaved, mode="w", fileformat="fastq")

    for read1, read2 in reader:
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

        # Progress reporting
        progress.update()

    reader.close()
    writer.close()
    logger.info(f'Finished')


def add_arguments(parser):
    # Files
    parser.add_argument("input", nargs='*',
                        help="Input paired fastq for read1 and read2. Use '-' if using stdin.")
    parser.add_argument("-o1", "--output1", default=None,
                        help="Output file name for read1. If not specified the result is written to stdout"
                             " as interleaved fastq.")
    parser.add_argument("-o2", "--output2", default=None,
                        help="Output file name for read2. If not specified but -o1 given the result is written"
                             " as interleaved fastq to o1.")
    # Options
    parser.add_argument("--interleaved", default=False, action='store_true',
                        help="Input is interleaved. ")
