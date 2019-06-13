"""
Take trimmed read barcodes sequences from headers (@HEADER_bc-seq)
and write FASTA files with unique barcodes
"""
import sys

import blr.utils as BLR


def main(args):
    """Takes a fastq file barcode sequences in the header and writes a barcode fasta file with only unique entries. """

    # Check python3 is being run
    if not BLR.pythonVersion(args.force_run): sys.exit()

    #
    # Filtering
    #
    if not args.filter == 1: BLR.report_progress('Filtering barcodes with less than ' + str(args.filter) + ' reads\n')

    #
    # Data processing
    #

    # Reading file and building initial bc dict with read counts
    bc_dict = dict()
    generator = BLR.FileReader(args.input_fastq)
    for read in generator.fastqReader():
        bc = read.header.split()[0].split('_')[-1]
        if bc in bc_dict:
            bc_dict[bc] += 1
        else:
            bc_dict[bc] = 1
    generator.close()

    # Indexing mode output writing
    bc_written = int()
    if args.index:
        # Make directory to put indexing files in
        index_dict, not_ATGC_index = reduceComplexity(bc_dict, args.index)
        try:
            import os
            os.mkdir(args.output_fasta)
        except FileExistsError:
            pass

        # Write one file per index
        for index in index_dict.keys():
            bc_id = int()
            output = args.output_fasta + '/' + str(index) + '.fa'
            with open(output, 'w') as openout:
                for barcode, read_count in index_dict[index].items():
                    if read_count < args.filter: continue
                    bc_id += 1
                    openout.write('>' + str(bc_id) + ':' + str(read_count) + ':' + str(barcode) + '\n' + str(barcode) + '\n')
                    bc_written += 1

    # Non-indexing mode output writing
    else:
        bc_id = int()
        output = args.output_fasta
        with open(output, 'w') as openout:
            for barcode, read_count in bc_dict.items():
                if read_count < args.filter: continue
                bc_id += 1
                openout.write('>' + str(bc_id) + ':' + str(read_count) + ':' + str(barcode) + '\n' + str(barcode) + '\n')
                bc_written += 1

    # Reporting
    BLR.report_progress('Unique BC count in input\t' + str(len(bc_dict.keys())))
    BLR.report_progress('Unique BC count in output\t' + str(bc_written))
    if args.index: BLR.report_progress('BC count where N was in index (Omitted from tot. BC count):\t' + str(not_ATGC_index))
    BLR.report_progress('Finished')

def reduceComplexity(bc_dict, index):
    """ Uses r first bases as indexes and divides files accordingly."""

    # Generate dict with possibilities indexes
    bases_list = ['A','T','C','G']
    index_dict = {'':dict()}
    not_ATGC_index = int()

    # Repeat for lenth of i
    for i in range(index):
        # Extend every key...
        for key in index_dict.copy().keys():
            # ... with one of every base
            for base in bases_list:
                index_dict[key+base] = dict()
            # Remove non-extended key
            del index_dict[key]

    # Classify reads into indexes
    for barcode, count in bc_dict.items():
        try: index_dict[barcode[:index]][barcode] = count
        except KeyError:
            not_ATGC_index += 1

    return index_dict, not_ATGC_index


def add_arguments(parser):
    parser.add_argument("input_fastq",
                        help="Read file with barcode sequences as last element of accession row, separated "
                             "by and underline. Example: '@ACCESSION_AGGTCGTCGATC'. Also handles "
                             "'@ACCESSSION_AGGTCGTCGATC MORE_ACCESSION'.")
    parser.add_argument("output_fasta", help="Output file name with unique barcode sequences.")

    parser.add_argument("-f", "--filter", type=int, default=1,
                        help="Filter file for minimum amount of read pairs. DEFAULT: 1")
    parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                       "Not recommended due to different function "
                                                                       "names in python 2 and 3.")
    parser.add_argument("-i", "--index", type=int, help="Divide BC sequences into descrete files due to their (-i) "
                                                        "first bases. DEFAULT: None")
    parser.add_argument("-s", "--space_separation", action="store_true", help='If BC is separated by <space> ( ) '
                                                                              'instead of <underline> (_)')
