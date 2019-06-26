"""
Removes barcode tags present at more than -M loci (corresponding to removing barcode tags from reads origin to droplets
which had more than -M molecules in one and the same droplet).
"""

import pysam
import logging

from tqdm import tqdm


logger = logging.getLogger(__name__)


def main(args):
    summary = Summary()
    allMolecules = AllMolecules()

    logger.info(f'Running analysis with {"{:,}".format(args.window)} bp window size')
    logger.info('Fetching reads')
    with pysam.AlignmentFile(args.x2_bam, 'rb') as infile:
        # Bam file stats
        prev_chrom = infile.references[0]
        summary.reads = infile.mapped + infile.unmapped
        for read in tqdm(infile.fetch(until_eof=True)):

            # Fetches barcode and genomic position. Position will be formatted so start < stop.
            BC_id, read_start, read_stop, summary = fetch_and_format(read, args.barcode_tag, summary=summary)
            if BC_id == None or read_start == 'unmapped': continue

            # Commit molecules between chromosomes
            if not prev_chrom == read.reference_name:
                allMolecules.reportAndRemoveAll(summary=summary)
                prev_chrom = read.reference_name

            if BC_id in allMolecules.cache_dict:
                molecule = allMolecules.cache_dict[BC_id]

                # Read is within args.window => add read to molecule.
                if (molecule.stop+args.window) >= read_start and molecule.stop < read_start:
                    molecule.addRead(stop=read_stop, read_header=read.query_name)
                    allMolecules.cache_dict[BC_id] = molecule

                # Overlapping reads => If not overlapping to it's mate, discard read.
                elif molecule.stop >= read_start:
                    if read.query_name in molecule.read_headers:
                        molecule.addRead(stop=read_stop, read_header=read.query_name)
                        allMolecules.cache_dict[BC_id] = molecule
                    else:
                        summary.overlapping_reads_in_pb += 1

                # Read is not within window => report old and initiate new molecule for that barcode.
                else:
                    allMolecules.report(molecule=molecule, summary=summary)
                    allMolecules.terminate(molecule=molecule)
                    molecule = Molecule(barcode=BC_id, start=read_start, stop=read_stop, read_header=read.query_name)
                    allMolecules.cache_dict[molecule.barcode] = molecule

            else:
                molecule = Molecule(barcode=BC_id, start=read_start, stop=read_stop, read_header=read.query_name)
                allMolecules.cache_dict[molecule.barcode] = molecule

    # Commit last chr molecules and log stats
    allMolecules.reportAndRemoveAll(summary=summary)
    summary.non_analyzed_reads = summary.unmapped_reads + summary.non_tagged_reads + summary.overlapping_reads_in_pb
    logger.info('Molecules analyzed')

    # Stats to output files and stdout
    if args.print_stats:
        summary.writeResultFiles(output_prefix=args.print_stats, threshold=args.threshold, filter_bam=args.output, Max_molecules=args.Max_molecules)
        summary.printStats(barcode_tag=args.barcode_tag, threshold=args.threshold, filter_bam=args.output)

    # Writes output bam file if wanted
    with pysam.AlignmentFile(args.x2_bam, 'rb') as openin:

        # OPEN FILES
        openfiles = dict()
        # IF SPLITTING INTO SEVERAL OUTPUTS
        if args.split:
            openfiles['no_bc'] = pysam.AlignmentFile(args.output + '.no_bc.bam', 'wb', template=openin)
            openfiles['not_phased'] = pysam.AlignmentFile(args.output + '.not_phased.bam', 'wb', template=openin)
            for reads_per_mol in range(1, args.Max_molecules+1):
                openfiles[reads_per_mol] = pysam.AlignmentFile(
                    args.output + '.' + str(reads_per_mol) + '_rpm.bam', 'wb', template=openin)
        # NORMAL FILTERED FILE
        else:
            open_file_key = 'all_reads'
            openfiles[open_file_key] = pysam.AlignmentFile(args.output, 'wb', template=openin)
            openout = openfiles[open_file_key]

        for read in tqdm(openin.fetch(until_eof=True)):

            # If no bc_id, just write it to out
            try: BC_id = read.get_tag(args.barcode_tag)
            except KeyError:
                BC_id = False

                # IF SPLITTING INTO SEVERAL OUTPUTS
                if args.split:
                    openout = openfiles['no_bc']

            # If too many molecules in cluster, change tag and header of read
            if BC_id in summary.barcode_removal_set:
                tmp_header_list = read.query_name.split('_')
                read.query_name = str(tmp_header_list[0]) + '_' + str(tmp_header_list[1])
                read.set_tag(args.barcode_tag, 'FILTERED', value_type='Z')
                summary.reads_with_removed_barcode += 1

                # IF SPLITTING INTO SEVERAL OUTPUTS
                if args.split:
                    openout = openfiles['no_bc']

            # IF SPLITTING INTO SEVERAL OUTPUTS
            elif args.split and BC_id:
                if not BC_id in summary.bc_to_numberMolOverReadThreshold:
                    openout = openfiles['not_phased']
                else:
                    mol_per_barcode = summary.bc_to_numberMolOverReadThreshold[BC_id]
                    openout = openfiles[mol_per_barcode]

            openout.write(read)

        # CLOSE FILES
        for openout in openfiles.values():
            openout.close()

    try:
        logger.info(f'Reads with barcodes removed:\t{"{:,}".format(summary.reads_with_removed_barcode)}\t'
                    f'({"%.2f" % ((summary.reads_with_removed_barcode/summary.reads)*100)}%)')
    except ZeroDivisionError:
        logger.warning('No reads passing filters found in file.')

def fetch_and_format(read, barcode_tag, summary):
    """
    Fetches barcode tag and turns read positions so read start < read stop
    """

    try: BC_id = read.get_tag(barcode_tag)
    except KeyError:
        summary.non_tagged_reads += 1
        BC_id = None

    if not read.is_unmapped:
        pos = (read.get_reference_positions()[0], read.get_reference_positions()[-1])
        read_start = min(pos)
        read_stop = max(pos)
    else:
        read_start, read_stop = 'unmapped', 'unmapped'
        summary.unmapped_reads += 1

    return BC_id, read_start, read_stop, summary

class Molecule:
    """
    Splits reads for barcodes into several molecules based on mapping proximity. Equivalent to several molecules being
    barcoded simultaneously in the same emulsion droplet (meaning with the same barcode).
    """

    def __init__(self, barcode, start, stop, read_header):
        """
        Creates a molecule defined by one read.
        """

        self.barcode = barcode
        self.start = start
        self.stop = stop
        self.read_headers = set()
        self.read_headers.add(read_header)

        self.number_of_reads = 1

    def addRead(self, stop, read_header):
        """
        Updates molecule's stop poistion, number of reads and header name set()
        """

        self.stop = stop
        self.read_headers.add(read_header)

        self.number_of_reads += 1


class AllMolecules:
    """
    Tracks all molecule information, with finsished molecules in .final_dict and molecules which still might get more
    reads in .cache_dict.
    """

    def __init__(self):
        self.cache_dict = dict()
        self.final_dict = dict()

    def report(self, molecule, summary):
        """
        Commit molecule to a set inside of an overlying dict where barcode:set(molecule1, molecule2...)
        """

        if not molecule.barcode in self.final_dict:
            self.final_dict[molecule.barcode] = set()

        self.final_dict[molecule.barcode].add(molecule)
        summary.molecules += 1

        return summary

    def terminate(self, molecule):
        """
        Removes a specific molecule from .cache_dict
        """

        del self.cache_dict[molecule.barcode]

    def reportAndRemoveAll(self, summary):
        """
        Commit all .cache_dict molecules to .final_dict and empty .cache_dict.
        """

        for molecule in self.cache_dict.values():
            self.report(molecule=molecule, summary=summary)
        self.cache_dict = dict()

        return summary

class Summary:

    def __init__(self):

        # Stats
        self.reads = int()
        self.molecules_over_threshold = int()
        self.molecules_result_dict = dict()
        self.molecules = int()
        self.non_tagged_reads = int()
        self.drops_without_molecules_over_threshold = int()
        self.overlapping_reads_in_pb = int()
        self.barcode_removal_set = set()
        self.reads_with_removed_barcode = int()
        self.unmapped_reads = int()
        self.bp_btw_reads = dict()
        self.non_analyzed_reads = int()

        # Filter bam file
        self.mol_rmvd_outbam = int()
        self.bc_rmvd_outbam = int()

        # Stats tracker needed to split bam files into separate according barcode per molecule
        self.bc_to_numberMolOverReadThreshold = dict()

    def printStats(self, barcode_tag, threshold, filter_bam):

        # Read stats
        logger.info('- Read stats -')
        logger.info(f'Total Reads in file:\t{"{:,}".format(self.reads)}')
        logger.info('- Reads skipped in analysis -')
        logger.info(f'Unmapped:\t{"{:,}".format(self.unmapped_reads)}')
        logger.info(f'Without {barcode_tag} tag:\t{"{:,}".format(self.non_tagged_reads)}')
        logger.info(f'Overlapping with other reads in molecule:\t{"{:,}".format(self.overlapping_reads_in_pb)}')
        logger.info('- Remaining reads -')
        logger.info(f'Reads analyzed:\t{"{:,}".format(self.reads-self.non_analyzed_reads)}')

        # Molecule stats
        logger.info('- Molecule stats -')
        logger.info(f'Molecules total:\t{"{:,}".format(self.molecules)}')
        logger.info(f'Molecules kept for stats (min read: {threshold}):\t{"{:,}".format(self.molecules_over_threshold)}')
        logger.info(f'BC consequently removed:\t{"{:,}".format(self.drops_without_molecules_over_threshold)}')

        # Filtering stats
        if filter_bam:
            logger.info('- Bam output stats -')
            logger.info(f'Molecules removed in output:\t{"{:,}".format(self.mol_rmvd_outbam)}')
            logger.info(f'BC removed in output:\t{"{:,}".format(self.bc_rmvd_outbam)}')

    def writeResultFiles(self, output_prefix, threshold, filter_bam, Max_molecules):

        # Opening all files
        molecules_per_bc_out = open((output_prefix + '.molecules_per_bc'), 'w')
        percent_bases_read = open((output_prefix + '.percent_bases_read'), 'w')
        reads_per_molecule_out = open((output_prefix + '.reads_per_molecule'), 'w')
        molecule_len_out = open((output_prefix + '.molecule_lengths'), 'w')
        everything = open((output_prefix + '.everything'), 'w')

        # Writing molecule-dependant stats
        for barcode_id in tqdm(self.molecules_result_dict.keys()):

            molecules_in_cluster = 0
            everything_cache_row = list()
            for molecule in self.molecules_result_dict[barcode_id]:

                reads_per_molecule_out.write(str(molecule[3]) + '\n')

                # Filter molecules with less than N reads (--threshold)
                if molecule[3] >= threshold:
                    molecules_in_cluster += 1
                    percent_bases_read.write(str(molecule[4]))
                    self.molecules_over_threshold += 1
                    molecule_len_out.write(str(molecule[2]) + '\n')
                    everything_cache_row.append((str(molecule[3]) + '\t' + str(molecule[2]) + '\t' + str(molecule[4]) + '\t' + str(molecule[4]) + '\t'  + str(barcode_id) + '\t'))

            # Skips clusters which as a result of -t does not have any molecules left.
            if molecules_in_cluster == 0:
                self.drops_without_molecules_over_threshold += 1
            else:

                # Stats tracker needed to split bam files into separate according barcode per molecule
                if not barcode_id in self.bc_to_numberMolOverReadThreshold: self.bc_to_numberMolOverReadThreshold[
                    barcode_id] = molecules_in_cluster

                molecules_per_bc_out.write(str(molecules_in_cluster) + '\t')
                if filter_bam:
                    if molecules_in_cluster > Max_molecules:
                        self.bc_rmvd_outbam += 1
                        self.mol_rmvd_outbam += molecules_in_cluster
                        self.barcode_removal_set.add(barcode_id)

                # Writes everything file afterwards in chunks (since it needs molecule per droplet)
                for row in everything_cache_row:
                    everything.write(row + '\t' + str(molecules_in_cluster) + '\n')

        # Close files
        for output_file in (molecules_per_bc_out, percent_bases_read, reads_per_molecule_out, molecule_len_out):
            output_file.close()

        # Distance between reads
        with open(output_prefix + '.lengths_between_readpairs', 'w') as openout:
            for length, number_of_times in self.bp_btw_reads.items():
                for i in range(number_of_times):
                    openout.write(str(length) + '\n')


def add_arguments(parser):
    parser.add_argument("x2_bam", help=".bam file tagged with BC:Z:<int> tags. Needs to be indexed, sorted & have "
                                       "duplicates removed.")
    parser.add_argument("output", help="Output filtered file.")
    parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                       "Not recommended due to different function "
                                                                       "names in python 2 and 3. DEFAULT: False")
    parser.add_argument("-t","--threshold", metavar='<INTEGER>', type=int, default=4,
                        help="Threshold for how many reads are required for including given molecule in statistics "
                             "(except_reads_per_molecule). DEFAULT: 4")
    parser.add_argument("-w", "--window", metavar='<INTEGER>', type=int, default=30000,
                        help="Window size cutoff for maximum distance "
                                                                              "in between two reads in one molecule. "
                                                                              "DEFAULT: 30000")
    parser.add_argument("-bc", "--barcode_tag", metavar='<STRING>', type=str, default='BC',
                        help="Bam file tag where barcode is stored. DEFAULT: BC")
    parser.add_argument("-ps", "--print_stats", metavar='<PREFIX>', type=str,
                        help="Write barcode statistics files (reads per molecule, molecule per barcode, molecule "
                             "lengths & coupling effiency). DEFAULT: None")
    parser.add_argument("-M", "--Max_molecules", metavar='<INTEGER>', type=int, default=500,
                        help="When using -f (--filter) this will remove barcode tags for those clusters which have more "
                             "than -M molecules. DEFAULT: 500")
    parser.add_argument("-s", "--split", action="store_true", help="Will intead of writing one output filtered file "
                                                                   "split output into separate files based on "
                                                                   "#mol/bc. The -f (--filter) output file name will instead be "
                                                                   "used as prefix for the output. -f option "
                                                                   "required. DEFAULT: None")
