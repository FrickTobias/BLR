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

    # Build allMolecules.final_dict[barcode][moleculeID] = molecule
    # molecule are instances of Molecule, defined as proximal reads if they are more than min_reads argument.
    logger.info(f'Running analysis with {"{:,}".format(args.window)} bp window size')
    logger.info('Fetching reads')
    with pysam.AlignmentFile(args.x2_bam, 'rb') as infile:
        allMolecules = build_molecule_dict(pysam_openfile=infile, barcode_tag=args.barcode_tag, window=args.window, min_reads=args.threshold, summary=summary)
        allMolecules.reportAndRemoveAll()

        summary.tot_reads = infile.mapped + infile.unmapped
        summary.unmapped_reads = infile.unmapped
        summary.mapped_reads = infile.mapped

    # Commit last chr molecules and log stats
    summary.non_analyzed_reads = summary.unmapped_reads + summary.non_tagged_reads + summary.overlapping_reads_in_pb
    logger.info('Molecules analyzed')

    # Writes output bam file if wanted
    with pysam.AlignmentFile(args.x2_bam, 'rb') as openin:
        with pysam.AlignmentFile(args.output, 'wb', template=openin) as openout:
            for read in tqdm(openin.fetch(until_eof=True)):

                # If no bc_id, just write it to out
                try: BC_id = read.get_tag(args.barcode_tag)
                except KeyError:
                    BC_id = False

                # If BC_id is not in allMolecules there the barcode does not have enough proximal reads to make a single molecule
                if BC_id in allMolecules.final_dict:

                    # If too many molecules in cluster, change tag and header of read
                    if len(allMolecules.final_dict[BC_id]) > args.Max_molecules:
                        tmp_header_list = read.query_name.split('_')
                        read.query_name = str(tmp_header_list[0]) + '_' + str(tmp_header_list[1])
                        read.set_tag(args.barcode_tag, 'FILTERED', value_type='Z')
                        summary.reads_with_removed_barcode += 1
                        if not BC_id in summary.barcode_removal_set:
                            summary.barcode_removal_set.add(BC_id)
                            summary.number_removed_molecules += len(allMolecules.final_dict[BC_id])

                openout.write(read)

    # Stats to output files and stdout
    summary.printStats(barcode_tag=args.barcode_tag, threshold=args.threshold, allMolecules=allMolecules)
    if args.print_stats:
        summary.writeMoleculeStats(output_prefix=args.print_stats, Max_molecules=args.Max_molecules,
                                 allMolecules=allMolecules)

def build_molecule_dict(pysam_openfile, barcode_tag, window, min_reads, summary):

    allMolecules = AllMolecules(min_reads=min_reads)

    prev_chrom = pysam_openfile.references[0]
    for read in tqdm(pysam_openfile.fetch(until_eof=True)):

        # Fetches barcode and genomic position. Position will be formatted so start < stop.
        BC_id, read_start, read_stop, summary = fetch_and_format(read, barcode_tag, summary=summary)
        if BC_id == None or read.is_unmapped: continue

        # Commit molecules between chromosomes
        if not prev_chrom == read.reference_name:
            allMolecules.reportAndRemoveAll()
            prev_chrom = read.reference_name

        if BC_id in allMolecules.cache_dict:
            molecule = allMolecules.cache_dict[BC_id]

            # Read is within window => add read to molecule (don't include overlapping reads).
            if (molecule.stop + window) >= read_start:
                if molecule.stop >= read_start and not read.query_name in molecule.read_headers:
                    summary.overlapping_reads_in_pb += 1
                else:
                    molecule.addRead(stop=read_stop, read_header=read.query_name)
                    allMolecules.cache_dict[BC_id] = molecule

            # Read is not within window => report old and initiate new molecule for that barcode.
            else:
                allMolecules.report(molecule=molecule)
                allMolecules.terminate(molecule=molecule)
                molecule = Molecule(barcode=BC_id, start=read_start, stop=read_stop, read_header=read.query_name)
                allMolecules.cache_dict[molecule.barcode] = molecule

        else:
            molecule = Molecule(barcode=BC_id, start=read_start, stop=read_stop, read_header=read.query_name)
            allMolecules.cache_dict[molecule.barcode] = molecule

    return allMolecules

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
        read_start = None
        read_stop = None

    return BC_id, read_start, read_stop, summary

class Molecule:
    """
    Splits reads for barcodes into several molecules based on mapping proximity. Equivalent to several molecules being
    barcoded simultaneously in the same emulsion droplet (meaning with the same barcode).
    """

    molecule_counter = int()

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

        Molecule.molecule_counter += 1
        self.ID = Molecule.molecule_counter

    def length(self):
        return(self.stop - self.start)

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

    def __init__(self, min_reads):

        # Min required reads for calling proximal reads a molecule
        self.min_reads = min_reads

        # Molecule tracking system
        self.cache_dict = dict()
        self.final_dict = dict()

    def report(self, molecule):
        """
        Commit molecule to a set inside of an overlying dict where barcode:set(molecule1, molecule2...)
        """

        if molecule.number_of_reads >= self.min_reads:
            if not molecule.barcode in self.final_dict:
                self.final_dict[molecule.barcode] = set()
            self.final_dict[molecule.barcode].add(molecule)

    def terminate(self, molecule):
        """
        Removes a specific molecule from .cache_dict
        """

        del self.cache_dict[molecule.barcode]

    def reportAndRemoveAll(self):
        """
        Commit all .cache_dict molecules to .final_dict and empty .cache_dict.
        """

        for molecule in self.cache_dict.values():
            self.report(molecule=molecule)
        self.cache_dict = dict()


class Summary:

    def __init__(self):

        self.mapped_reads = int()

        # Stats
        self.tot_reads = int()
        self.non_tagged_reads = int()
        self.overlapping_reads_in_pb = int()
        self.barcode_removal_set = set()
        self.reads_with_removed_barcode = int()
        self.unmapped_reads = int()
        self.non_analyzed_reads = int()
        self.number_removed_molecules = int()

        # Stats tracker needed to split bam files into separate according barcode per molecule
        self.bc_to_numberMolecules = dict()

    def printStats(self, barcode_tag, threshold, allMolecules):

        # Read stats
        logger.info('- Read stats -')
        logger.info(f'Total Reads in file:\t{"{:,}".format(self.tot_reads)}')
        logger.info('- Reads skipped in analysis -')
        logger.info(f'Unmapped:\t{"{:,}".format(self.unmapped_reads)}')
        logger.info(f'Without {barcode_tag} tag:\t{"{:,}".format(self.non_tagged_reads)}')
        logger.info(f'Overlapping with other reads in molecule:\t{"{:,}".format(self.overlapping_reads_in_pb)}')
        logger.info('- Remaining reads -')
        logger.info(f'Reads analyzed:\t{"{:,}".format(self.tot_reads - self.non_analyzed_reads)}')

        # Molecule stats
        logger.info('- Molecule stats -')
        logger.info(f'Molecules total (min read {threshold}):\t{"{:,}".format(sum(len(all) for all in allMolecules.final_dict.values()))}')
        logger.info(f'Barcodes removed:\t{len(self.barcode_removal_set)}')
        logger.info(f'Molecules removed:\t{self.number_removed_molecules}')

        try:
            logger.info(f'Reads with barcodes removed:\t{"{:,}".format(self.reads_with_removed_barcode)}\t'
                        f'({"%.2f" % ((summary.reads_with_removed_barcode/self.tot_reads)*100)}%)')
        except ZeroDivisionError:
            logger.warning('No reads passing filters found in file.')

    def writeMoleculeStats(self, output_prefix, allMolecules):

        # Opening all files
        molecules_per_bc_out = open((output_prefix + '.molecules_per_bc'), 'w')
        reads_per_molecule_out = open((output_prefix + '.reads_per_molecule'), 'w')
        molecule_len_out = open((output_prefix + '.molecule_lengths'), 'w')
        everything = open((output_prefix + '.everything'), 'w')

        # Writing molecule-dependant stats
        for barcode in tqdm(allMolecules.final_dict):
            molecules_per_bc_out.write(str(len(allMolecules.final_dict[barcode])))
            for molecule in (allMolecules.final_dict[barcode]):
                everything.write(str(molecule.number_of_reads) + '\t' + str(molecule.length()) + '\t' + str(barcode) + '\t' + str(len(allMolecules.final_dict[barcode])) + '\n')

            # Stats tracker needed to split bam files into separate according barcode per molecule
            self.bc_to_numberMolecules[barcode] = len(allMolecules.final_dict[barcode])

        # Close files
        for output_file in (molecules_per_bc_out, reads_per_molecule_out, molecule_len_out, everything):
            output_file.close()

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
