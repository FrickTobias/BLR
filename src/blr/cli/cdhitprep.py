"""
Take trimmed read barcodes sequences from headers (@HEADER_bc-seq)
and write FASTA files with unique barcodes
"""
import os
import logging

import dnaio
from collections import defaultdict
from itertools import product

logger = logging.getLogger(__name__)


def main(args):
    """
    Takes a fastq file barcode sequences in the header and writes a barcode
    fasta file with only unique entries.
    """

    logger.info(f'Filtering barcodes with less than {args.filter} reads')

    # Reading file and building initial bc dict with read counts
    barcode_counts = defaultdict(int)
    separator = "_" if not args.space_separation else " "
    with dnaio.open(args.input_fastq, fileformat="fastq", mode="r") as reader:
        for read in reader:
            barcode_sequence = read.name.split()[0].split(separator)[-1]
            barcode_counts[barcode_sequence] += 1

    # Indexing mode output writing
    if args.index:

        # Get barcode counts for each index of length=args.index.
        indexed_barcode_count, not_atcg_index = reduce_complexity(barcode_counts, index_size=args.index)

        # Make directory to put indexing files in unless already present
        try:
            os.mkdir(args.output_fasta)
        except FileExistsError:
            pass

        # Write one file per index
        for index_sequence in indexed_barcode_count.keys():
            output = f'{args.output_fasta}/{index_sequence}.fa'

            logger.info(f'Writing output to {output}')

            with dnaio.open(output, fileformat="fasta", mode='w') as openout:
                for bc_id, (barcode, read_count) in enumerate(indexed_barcode_count[index_sequence].items(), start=1):
                    if read_count < args.filter:
                        continue

                    fasta_name = f'>{bc_id}:{read_count}:{barcode}'
                    fasta_entry = dnaio.Sequence(name=fasta_name, sequence=barcode)
                    openout.write(fasta_entry)

    # Non-indexing mode output writing
    else:
        # Check if file format matches fasta
        if any(args.output_fasta.endswith(extension) for extension in ['.fa', '.fasta']):
            output = args.output_fasta
        else:
            output = f'{args.output_fasta}.fasta'

        logger.info(f'Writing output to {output}')

        # Write all output to one file.
        with dnaio.open(output, fileformat="fasta", mode="w") as openout:
            for bc_id, (barcode, read_count) in enumerate(barcode_counts.items(), start=1):
                if read_count < args.filter:
                    continue

                fasta_name = f'>{bc_id}:{read_count}:{barcode}'
                fasta_entry = dnaio.Sequence(name=fasta_name, sequence=barcode)
                openout.write(fasta_entry)

    # Reporting
    logger.info(f'Unique BC count in input:\t{len(barcode_counts)}')
    if args.index:
        logger.info(f'BC count where N was in index (Omitted from tot. BC count):\t{not_atcg_index}')
    logger.info("Finished")


def reduce_complexity(barcode_counts, index_size=1):
    """
    Uses first bases as indexes and divides files accordingly to reduce complexity.
    :param barcode_counts: Dictionary containing barcode sequences as key and counts as values
    :param index_size: Number of bases to use for complexity reduction. Reduction is proportional to 4^(index_size)
    """

    # Generate dict with possibilities indexes
    indeces = product('ATCG', repeat=index_size)
    indexed_barcode_count = {''.join(index): dict() for index in indeces}

    not_atcg_index = int()

    # Classify reads into indexes
    for barcode, count in barcode_counts.items():
        try:
            indexed_barcode_count[barcode[:index_size]][barcode] = count
        except KeyError:
            not_atcg_index += 1

    return indexed_barcode_count, not_atcg_index


def add_arguments(parser):
    parser.add_argument("input_fastq",
                        help="Read file with barcode sequences as last element of accession row, separated "
                             "by and underline. Example: '@ACCESSION_AGGTCGTCGATC'. Also handles "
                             "'@ACCESSSION_AGGTCGTCGATC MORE_ACCESSION'.")
    parser.add_argument("output_fasta", help="Output file name with unique barcode sequences.")

    parser.add_argument("-f", "--filter", type=int, default=1,
                        help="Filter file for minimum amount of read pairs. DEFAULT: 1")
    parser.add_argument("-i", "--index", type=int, default=None, help="Divide BC sequences into descrete files due "
                                                                      "to their (-i) first bases. DEFAULT: None")
    parser.add_argument("-s", "--space_separation", action="store_true", help='If BC is separated by <space> ( ) '
                                                                              'instead of <underline> (_)')
