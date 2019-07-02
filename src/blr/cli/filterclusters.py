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
    logger.info(f"Running analysis with {args.window:,} bp window size")

    # Build molecule dictionary used for counting #reads/molecule & #molecules/barcode for filtering of bam file
    with pysam.AlignmentFile(args.x2_bam, "rb") as infile:
        molecule_dict = build_molecule_dict(pysam_openfile=infile, barcode_tag=args.barcode_tag, window=args.window,
                                            min_reads=args.threshold, summary=summary)

        summary.tot_reads = infile.mapped + infile.unmapped
        summary.unmapped_reads = infile.unmapped
        summary.mapped_reads = infile.mapped
        summary.non_analyzed_reads = summary.unmapped_reads + summary.non_tagged_reads + summary.overlapping_reads_in_pb

    # Writes filtered out
    with pysam.AlignmentFile(args.x2_bam, "rb") as openin, \
            pysam.AlignmentFile(args.output, "wb", template=openin) as openout:
        logger.info("Writing filtered bam file")
        for read in tqdm(openin.fetch(until_eof=True)):
            barcode = fetch_bc(pysam_read=read, barcode_tag=args.barcode_tag)

            # If barcode is not in allMolecules the barcode does not have enough proximal reads to make a single
            # molecule. If the barcode has more than <max_molecules> molecules, remove it from the read.
            if barcode in molecule_dict and len(molecule_dict[barcode]) > args.max_molecules:
                read = strip_barcode(pysam_read=read ,barcode_tag=args.barcode_tag)

                summary.reads_with_removed_barcode += 1
                if not barcode in summary.removed_barcodes:
                    summary.removed_barcodes.add(barcode)
                    summary.removed_molecules += len(molecule_dict[barcode])

            openout.write(read)

    summary.print_stats(barcode_tag=args.barcode_tag, molecule_dict=molecule_dict)

    # Write molecule/barcode file stats
    if args.stats_file:
        logger.info("Writing statistics files")
        summary.write_molecule_stats(output_prefix=args.stats_file, molecule_dict=molecule_dict)

def build_molecule_dict(pysam_openfile, barcode_tag, window, min_reads, summary):
    """
    Build allMolecules.final_dict, [barcode][moleculeID] = molecule where molecule are instances of Molecule object,
    defined as more than <min_reads> within <window> from each other.
    :param pysam_openfile: Pysam open file instance.
    :param barcode_tag: Tag used to store barcode in bam file (usually BC).
    :param window: Max distance between reads to include in the same molecule.
    :param min_reads: Minimum reads to include molecule in allMolecules.final_dict
    :param summary: Custom summary intsance
    :return: dict[barcode][molecule] = moleculeInstance
    """

    allMolecules = AllMolecules(min_reads=min_reads)

    prev_chrom = pysam_openfile.references[0]
    logger.info("Dividing barcodes into molecules")
    for read in tqdm(pysam_openfile.fetch(until_eof=True)):

        # Fetches barcode and genomic position. Position will be formatted so start < stop.
        barcode = fetch_bc(pysam_read=read, barcode_tag=barcode_tag, summary=summary)
        if barcode and read.is_unmapped == False:
            read_start, read_stop = sorted((read.get_reference_positions()[0], read.get_reference_positions()[-1]))

            # Commit molecules between chromosomes
            if not prev_chrom == read.reference_name:
                allMolecules.report_and_remove_all()
                prev_chrom = read.reference_name

            if barcode in allMolecules.cache_dict:
                molecule = allMolecules.cache_dict[barcode]

                # Read is within window => add read to molecule (don't include overlapping reads).
                if (molecule.stop + window) >= read_start:
                    if molecule.stop >= read_start and not read.query_name in molecule.read_headers:
                        summary.overlapping_reads_in_pb += 1
                    else:
                        molecule.add_read(stop=read_stop, read_header=read.query_name)
                        allMolecules.cache_dict[barcode] = molecule

                # Read is not within window => report old and initiate new molecule for that barcode.
                else:
                    allMolecules.report(molecule=molecule)
                    allMolecules.terminate(molecule=molecule)
                    molecule = Molecule(barcode=barcode, start=read_start, stop=read_stop, read_header=read.query_name)
                    allMolecules.cache_dict[molecule.barcode] = molecule

            else:
                molecule = Molecule(barcode=barcode, start=read_start, stop=read_stop, read_header=read.query_name)
                allMolecules.cache_dict[molecule.barcode] = molecule

    allMolecules.report_and_remove_all()

    return allMolecules.final_dict

def fetch_bc(pysam_read, barcode_tag, summary=None):
    """
    Fetches barcode from a bam file tag, returns None if reads isn't tagged.
    """

    try: barcode = pysam_read.get_tag(barcode_tag)
    except KeyError:
        barcode = None

        if summary:
            summary.non_tagged_reads += 1

    return barcode

def strip_barcode(pysam_read, barcode_tag):
    """
    Strips an alignment from its barcode sequence. Keeps information in header but adds FILTERED prior to bc info.
    """

    # Modify read instance
    header, bc_info = pysam_read.query_name.split("_", maxsplit=1)
    pysam_read.query_name = header + "_" + "FILTERED-" + bc_info
    pysam_read.set_tag(barcode_tag, "FILTERED", value_type="Z")

    return pysam_read

class Molecule:
    """
    A Splitting of barcode read groups into several molecules based on mapping proximity. Equivalent to several
    molecules being barcoded simultaneously in the same emulsion droplet (meaning with the same barcode).
    """

    molecule_counter = int()

    def __init__(self, barcode, start, stop, read_header):
        """
        :param barcode: barcode ID
        :param start: min(read_mapping_positions)
        :param stop: max(read_mapping_positions)
        :param read_header: read ID
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

    def add_read(self, stop, read_header):
        """
        Updates molecule's stop position, number of reads and header name set()
        """

        self.stop = stop
        self.read_headers.add(read_header)

        self.number_of_reads += 1


class AllMolecules:
    """
    Tracks all molecule information, with finished molecules in .final_dict, and molecules which still might get more
    reads in .cache_dict.
    """

    def __init__(self, min_reads):
        """
        :param min_reads: Minimum reads required to add molecule to .final_dict from .cache_dict
        """

        # Min required reads for calling proximal reads a molecule
        self.min_reads = min_reads

        # Molecule tracking system
        self.cache_dict = dict()
        self.final_dict = dict()

    def report(self, molecule):
        """
        Commit molecule to .final_dict, if molecule.reads >= min_reads
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

    def report_and_remove_all(self):
        """
        Commit all .cache_dict molecules to .final_dict and empty .cache_dict (provided they meet criterias by report
        function).
        """

        for molecule in self.cache_dict.values():
            self.report(molecule=molecule)
        self.cache_dict = dict()

class Summary:
    """
    Gathers all stats generated during analysis
    """

    def __init__(self):

        self.tot_reads = int()
        self.mapped_reads = int()

        self.unmapped_reads = int()
        self.non_tagged_reads = int()
        self.overlapping_reads_in_pb = int()
        self.non_analyzed_reads = int() # Sum of the three int() above

        self.removed_barcodes = set()
        self.removed_molecules = int()
        self.reads_with_removed_barcode = int()

    def print_stats(self, barcode_tag, molecule_dict):
        """
        Prints stats to terminal
        """

        # Read stats
        logger.info("- Read stats -")
        logger.info(f"Total Reads in file:\t{self.tot_reads:,}")
        logger.info("- Reads skipped in analysis -")
        logger.info(f"Unmapped:\t{self.unmapped_reads:,}")
        logger.info(f"Without {barcode_tag} tag:\t{self.non_tagged_reads:,}")
        logger.info(f"Overlapping with other reads in molecule:\t{self.overlapping_reads_in_pb:,}")
        logger.info("- Remaining reads -")
        logger.info(f"Reads analyzed:\t{self.tot_reads - self.non_analyzed_reads:,}")

        # Molecule stats
        logger.info("- Molecule stats -")
        logger.info(f"Molecules total:\t{sum(len(all) for all in molecule_dict.values()):,}")
        logger.info(f"Barcodes removed:\t{len(self.removed_barcodes)}")
        logger.info(f"Molecules removed:\t{self.removed_molecules}")

        try:
            logger.info(f"Reads with barcodes removed:\t{self.reads_with_removed_barcode:,}"
                        f"({self.reads_with_removed_barcode/self.tot_reads:.2%})")
        except ZeroDivisionError:
            logger.warning("No reads passing filters found in file.")

    def write_molecule_stats(self, output_prefix, molecule_dict):
        """
        Writes stats file for molecules and barcode with information like how many reads, barcodes, molecules etc they
        have
        """

        # Opening all files
        molecules_per_bc = open((output_prefix + ".molecules_per_bc"), "w")
        molecule_stats = open((output_prefix + ".molecule_stats"), "w")

        # Writing molecule-dependant stats
        for barcode in tqdm(molecule_dict):
            number_of_molecules = len(molecule_dict[barcode])
            molecules_per_bc.write(str(number_of_molecules))
            for molecule in (molecule_dict[barcode]):
                molecule_stats.write(str(molecule.number_of_reads) + "\t" + str(molecule.length()) + "\t" + barcode +
                                     "\t" + str(number_of_molecules) + "\n")

        # Close files
        for output_file in (molecules_per_bc, molecule_stats):
            output_file.close()

def add_arguments(parser):
    parser.add_argument("x2_bam", help=".bam file tagged with BC:Z:<int> tags. Needs to be indexed, sorted & have "
                                       "duplicates removed.")
    parser.add_argument("output", help="Output filtered file.")

    parser.add_argument("-t","--threshold", metavar="<INTEGER>", type=int, default=4,
                        help="Threshold for how many reads are required for including given molecule in statistics "
                             "(except_reads_per_molecule). DEFAULT: 4")
    parser.add_argument("-w", "--window", metavar="<INTEGER>", type=int, default=30000,
                        help="Window size cutoff for maximum distance in between two reads in one molecule. DEFAULT: "
                             "30000")
    parser.add_argument("-bc", "--barcode_tag", metavar="<STRING>", type=str, default="BC",
                        help="Bam file tag where barcode is stored. DEFAULT: BC")
    parser.add_argument("-s", "--stats_file", metavar="<PREFIX>", type=str,
                        help="Write barcode/molecule statistics files. DEFAULT: None")
    parser.add_argument("-M", "--max_molecules", metavar="<INTEGER>", type=int, default=500,
                        help="When using -f (--filter) this will remove barcode tags for those clusters which have more "
                             "than -M molecules. DEFAULT: 500")
