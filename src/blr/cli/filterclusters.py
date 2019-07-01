"""
Removes barcode tags present at more than -M loci (corresponding to removing barcode tags from reads origin to droplets
which had more than -M molecules).
"""

import pysam
import logging

import blr.utils as BLR

logger = logging.getLogger(__name__)


def main(args):
    summary = Summary()
    molecules = Molecules()
    prev_chrom = 'chr1'

    # Open file, loop over all reads
    logger.info(f'Running analysis with {"{:,}".format(args.window)} bp window size')
    logger.info('Fetching reads')
    progress = BLR.ProgressReporter('Reads processed', 1000000)
    with pysam.AlignmentFile(args.x2_bam, 'rb') as infile:
        for read in infile.fetch(until_eof=True):

            # Progress reporting. Upmost because of continue-statement in loop
            progress.update()

            # Fetches barcode and genomic position. Position will be formatted so start < stop.
            BC_id, read_start, read_stop, summary = fetch_and_format(read, args.barcode_tag, summary=summary)
            # If read is unmapped or does not have barcode, skip
            if BC_id == None or read_start == 'unmapped': continue

            # Commit molecules between chromosomes
            if not prev_chrom == read.reference_name:
                molecules.reportAndRemoveAll(summary=summary)
                prev_chrom = read.reference_name

            # If BC_id already has seen prior reads
            if BC_id in molecules.dictionary:
                molecule_stop = molecules.dictionary[BC_id]['stop']

                # Read is within window => add read to molecule.
                if (molecule_stop+args.window) >= read_start and molecule_stop < read_start:
                    molecules.addRead(name=BC_id, read_start=read_start, read_stop=read_stop, read_name=read.query_name, summary=summary)

                # Overlapping reads => If not overlapping to it's mate, discard read.
                elif molecule_stop >= read_start:
                    if read.query_name in molecules.dictionary[BC_id]['reads']:
                        molecules.addRead(name=BC_id, read_start=read_start, read_stop=read_stop, read_name=read.query_name, summary=summary)
                    else:
                        summary.overlapping_reads_in_pb += 1

                # Read is not within window => report old and initate new molecule.
                else:
                    summary.reportMolecule(name=BC_id, molecule=molecules.dictionary[BC_id])
                    molecules.terminate(name=BC_id)
                    molecules.initiate(name=BC_id, start=read_start, stop=read_stop, read_name=read.query_name, summary=summary)

            # No previous reads for this bc has been discovered
            else:
                molecules.initiate(name=BC_id, start=read_start, stop=read_stop, read_name=read.query_name, summary=summary)

    # Commit last chr molecules and log stats
    molecules.reportAndRemoveAll(summary=summary)
    summary.reads = progress.position
    summary.non_analyzed_reads = summary.unmapped_reads + summary.non_tagged_reads + summary.overlapping_reads_in_pb
    logger.info('Molecules analyzed')

    # Stats to output files and stdout
    summary.writeResultFiles(output_prefix=args.output_prefix, threshold=args.threshold, filter_bam=args.filter_bam, Max_molecules=args.Max_molecules)
    summary.printStats(barcode_tag=args.barcode_tag, threshold=args.threshold, filter_bam=args.filter_bam)

    # Writes output bam file if wanted
    if args.filter_bam:

        progressBar = BLR.ProgressBar(name='Writing filtered bam file', min=0, max=summary.reads, step=1)
        with pysam.AlignmentFile(args.x2_bam, 'rb') as openin:

            # OPEN FILES
            openfiles = dict()
            # IF SPLITTING INTO SEVERAL OUTPUTS
            if args.split:
                openfiles['no_bc'] = pysam.AlignmentFile(args.filter_bam + '.no_bc.bam', 'wb', template=openin)
                openfiles['not_phased'] = pysam.AlignmentFile(args.filter_bam + '.not_phased.bam', 'wb', template=openin)
                for reads_per_mol in range(1, args.Max_molecules+1):
                    openfiles[reads_per_mol] = pysam.AlignmentFile(
                        args.filter_bam + '.' + str(reads_per_mol) + '_rpm.bam', 'wb', template=openin)
            # NORMAL FILTERED FILE
            else:
                open_file_key = 'all_reads'
                openfiles[open_file_key] = pysam.AlignmentFile(args.filter_bam, 'wb', template=openin)
                openout = openfiles[open_file_key]

            for read in openin.fetch(until_eof=True):

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

                # Writ to out & update progress bar
                openout.write(read)
                progressBar.update()

            # CLOSE FILES
            for openout in openfiles.values():
                openout.close()

        # End progress bar
        progressBar.terminate()

    try:
        logger.info(f'Reads with barcodes removed:\t{"{:,}".format(summary.reads_with_removed_barcode)}\t'
                    f'({"%.2f" % ((summary.reads_with_removed_barcode/summary.reads)*100)}%)')
    except ZeroDivisionError:
        logger.warning('No reads passing filters found in file.')

def fetch_and_format(read, barcode_tag, summary):
    """

    :param read:
    :return:
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

class Molecules:
    """
    Tmp storage for molecules which might still get more reads assigned to them. Basically a dict with customised functions.
    """

    def __init__(self):
        self.dictionary = dict()

    def initiate(self, name, start, stop, read_name, summary):

        summary.molecules += 1

        self.dictionary[name] = dict()
        self.dictionary[name]['start'] = start
        self.dictionary[name]['stop'] = stop
        self.dictionary[name]['number_of_reads'] = 1
        self.dictionary[name]['bases_btw_inserts'] = 0
        self.dictionary[name]['bases_read'] = stop - start
        self.dictionary[name]['reads'] = set()
        self.dictionary[name]['reads'].add(read_name)

        return summary

    def addRead(self, name, read_start, read_stop, read_name, summary):

        # Tracks distances between read pairs, won't add value if it is mate to read
        if not read_name in self.dictionary[name]['reads']:
            bp_btw_reads = read_start - self.dictionary[name]['stop']

            # Tracking the different lengths
            if not bp_btw_reads in summary.bp_btw_reads:
                summary.bp_btw_reads[bp_btw_reads] = int()
            summary.bp_btw_reads[bp_btw_reads] += 1

        # Builds new object
        self.dictionary[name]['bases_read'] += read_stop - read_start
        self.dictionary[name]['bases_btw_inserts'] += read_start - self.dictionary[name]['stop']
        self.dictionary[name]['stop'] = read_stop
        self.dictionary[name]['number_of_reads'] += 1
        self.dictionary[name]['reads'].add(read_name)

        return summary

    def terminate(self, name):

        del self.dictionary[name]

    def reportAndRemoveAll(self, summary):

        for BC_id in self.dictionary.copy().keys():
            summary.reportMolecule(name=BC_id, molecule=self.dictionary[BC_id])
            del self.dictionary[BC_id]

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

    def reportMolecule(self, name, molecule):

        # Fetching and formatting information
        start = molecule['start']
        stop = molecule['stop']
        length = stop-start
        num_reads = molecule['number_of_reads']
        percent_bases_read = molecule['bases_read']/(molecule['bases_btw_inserts'] + molecule['bases_read'])

        # Tries to append to list of tuples, otherwise creates a tuple list as value for given barcode id
        try: self.molecules_result_dict[name]
        except KeyError:
            self.molecules_result_dict[name] = set()

        # Save in summary dictionary
        self.molecules_result_dict[name].add((start, stop, length, num_reads, percent_bases_read))

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
        progressBar = BLR.ProgressBar(name='Writing stats files', min=0, max = self.molecules, step = 1)

        # Writing molecule-dependant stats
        for barcode_id in self.molecules_result_dict.keys():

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

                progressBar.update()

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
                        self.mary.mol_rmvd_outbam += molecules_in_cluster
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

        progressBar.terminate()

def add_arguments(parser):
    parser.add_argument("x2_bam", help=".bam file tagged with @RG tags and duplicates removed. Needs to be indexed.")
    parser.add_argument("output_prefix", help="prefix for results files (.coupling_efficiency, "
                                              ".reads_per_molecule, .molecule_per_barcode and .molecule_lengths). "
                                              "Will also write a .everything file containing rows (equivalent to "
                                              "molecules) with all the information from the other files.")
    parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                       "Not recommended due to different function "
                                                                       "names in python 2 and 3. DEFAULT: False")
    parser.add_argument("-t","--threshold", metavar='<INTEGER>', type=int, default=4, help="Threshold for how many reads are required for including given molecule in statistics (except_reads_per_molecule). DEFAULT: 4")
    parser.add_argument("-w", "--window", metavar='<INTEGER>', type=int, default=30000, help="Window size cutoff for maximum distance "
                                                                              "in between two reads in one molecule. "
                                                                              "DEFAULT: 30000")
    parser.add_argument("-bc", "--barcode_tag", metavar='<STRING>', type=str, default='BC', help="Bam file tag where barcode is stored. DEFAULT: BC")
    parser.add_argument("-f", "--filter_bam", metavar='<OUTPUT>', type=str, default=False, help="Write an output bam file where reads have their barcode tag removed if it is "
                                                                            "from a cluster with more molecules than -M (--Max_molecules) molecules. DEFAULT: False")
    parser.add_argument("-M", "--Max_molecules", metavar='<INTEGER>', type=int, default=500, help="When using -f (--filter) this will remove barcode tags for those clusters which have more than -M molecules. DEFAULT: 500")
    parser.add_argument("-s", "--split", action="store_true", help="Will intead of writing one output filtered file "
                                                                   "split output into separate files based on "
                                                                   "#mol/bc. The -f (--filter) output file name will instead be "
                                                                   "used as prefix for the output. -f option "
                                                                   "required. DEFAULT: None")
