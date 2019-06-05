#! /usr/bin python3


def main():
    import blr.utils as BLR, sys

    #
    # Argument parsing
    #
    args = readArgs().parse()

    # Check python3 is being run
    if not BLR.pythonVersion(args.force_run): sys.exit()

    #
    # Process data
    #

    BLR.report_progress('Starting')
    progress = BLR.ProgressReporter('Read pairs processed', 1000000)
    generator = BLR.FileReader(args.r1, args.r2)
    with open(args.out_r1, 'w') as openr1, open(args.out_r2, 'w') as openr2:
        for read1, read2 in generator.fastqPairedReader():

            # Adjusting for BC
            bc_seq = read1.seq[:20]
            read1.seq = read1.seq[20:]
            read1.qual = read1.qual[20:]

            # Header parsing
            name_and_pos_r1 = read1.header.split()[0]
            read_and_index_r1 = read1.header.split()[1]
            name_and_pos_r2 = read2.header.split()[0]
            read_and_index_r2 = read2.header.split()[1]

            # Save header to read instances
            read1.header = name_and_pos_r1 + '_' + bc_seq + ' ' + read_and_index_r1
            read2.header = name_and_pos_r2 + '_' + bc_seq + ' ' + read_and_index_r2

            # Write to out
            openr1.write(read1.fastq_string())
            openr2.write(read2.fastq_string())

            # Progress reporting
            progress.update()

    generator.close()
    BLR.report_progress('Finished')


class readArgs:
    """
    Reads arguments and handles basic error handling like python version control etc.
    """

    def __init__(self):

        args = readArgs.parse(self)
        readArgs.pythonVersion(self, args)

    def parse(self):

        #
        # Imports & globals
        #
        import argparse

        parser = argparse.ArgumentParser(description="Extracts barcode sequences by moving 20 bp from 5' end "
                                                     "of the read to the header, separated by an underline.")

        # Arguments
        parser.add_argument("r1", help="Read 1 fastq file")
        parser.add_argument("r2", help="Read 2 fastq file")
        parser.add_argument("out_r1", help="Read 1 output")
        parser.add_argument("out_r2", help="Read 2 output")

        # Options
        parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                           "Not recommended due to different function "
                                                                           "names in python 2 and 3.")

        return parser.parse_args()

    def pythonVersion(self, args):
        """ Makes sure the user is running python 3."""

        #
        # Version control
        #
        import sys
        if sys.version_info.major == 3:
            pass
        else:
            sys.stderr.write('\nWARNING: you are running python ' + str(
                sys.version_info.major) + ', this script is written for python 3.')
            if not args.force_run:
                sys.stderr.write('\nAborting analysis. Use -F (--Force) to run anyway.\n')
                sys.exit()
            else:
                sys.stderr.write('\nForcing run. This might yield inaccurate results.\n')


if __name__=="__main__": main()
