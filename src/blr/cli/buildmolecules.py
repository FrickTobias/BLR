"""
Tags .bam file with molecule information based on barcode sequence and genomic proximity.

A molecule is defined by having 1) minimum --threshold reads and including all reads with the same barcode which are 2)
a maximum distance of --window between any given reads.
"""

import pysam
import logging

from tqdm import tqdm

logger = logging.getLogger(__name__)


def main(args):
    summary = Summary()

    # Build molecules from BCs and reads
    with pysam.AlignmentFile(args.bam, "rb") as infile:
        bc_to_mol_dict, header_to_mol_dict = build_molecules(pysam_openfile=infile,
                                                             barcode_tag=args.barcode_tag,
                                                             window=args.window, min_reads=args.threshold,
                                                             summary=summary)

    # Writes filtered out
    with pysam.AlignmentFile(args.bam, "rb") as openin, \
            pysam.AlignmentFile(args.output, "wb", template=openin) as openout:
        logger.info("Writing filtered bam file")
        for read in tqdm(openin.fetch(until_eof=True)):
            name = read.query_name

            # If barcode is not in all_molecules the barcode does not have enough proximal reads to make a single
            # molecule. If the barcode has more than <max_molecules> molecules, remove it from the read.
            if name in header_to_mol_dict:
                summary.reads_tagged += 1

                # Set tags
                molecule_ID = header_to_mol_dict[name]
                read.set_tag(args.molecule_tag, molecule_ID)
                bc_num_molecules = len(bc_to_mol_dict[read.get_tag(args.barcode_tag)])
                read.set_tag(args.number_tag, bc_num_molecules)

            openout.write(read)

    summary.print_stats()

    # Write molecule/barcode file stats
    if args.stats_file:
        logger.info("Writing statistics files")
        summary.write_molecule_stats(output_prefix=args.stats_file, molecule_dict=bc_to_mol_dict)


def build_molecules(pysam_openfile, barcode_tag, window, min_reads, summary):
    """
    Builds all_molecules.bc_to_mol ([barcode][moleculeID] = molecule) and
    all_molecules.header_to_mol ([read_name]=mol_ID)
    :param pysam_openfile: Pysam open file instance.
    :param barcode_tag: Tag used to store barcode in bam file (usually BC).
    :param window: Max distance between reads to include in the same molecule.
    :param min_reads: Minimum reads to include molecule in all_molecules.bc_to_mol
    :param summary: Custom summary intsance
    :return: dict[barcode][molecule] = moleculeInstance, dict[read_name] = mol_ID
    """

    all_molecules = AllMolecules(min_reads=min_reads)

    prev_chrom = pysam_openfile.references[0]
    logger.info("Dividing barcodes into molecules")
    for read in tqdm(pysam_openfile.fetch(until_eof=True)):
        summary.tot_reads += 1
        if read.is_duplicate:
            summary.duplicates += 1
            continue

        # Fetches barcode and genomic position. Position will be formatted so start < stop.
        barcode = fetch_bc(pysam_read=read, barcode_tag=barcode_tag, summary=summary)
        if barcode and not read.is_unmapped:
            read_start, read_stop = sorted((read.reference_start, read.reference_end))

            # Commit molecules between chromosomes
            if not prev_chrom == read.reference_name:
                all_molecules.report_and_remove_all()
                prev_chrom = read.reference_name

            if barcode in all_molecules.cache_dict:
                molecule = all_molecules.cache_dict[barcode]

                # Read is within window => add read to molecule (don't include overlapping reads).
                if (molecule.stop + window) >= read_start:
                    if molecule.stop >= read_start and read.query_name not in molecule.read_headers:
                        summary.overlapping_reads_in_molecule += 1
                    else:
                        molecule.add_read(stop=read_stop, read_header=read.query_name)
                        all_molecules.cache_dict[barcode] = molecule

                # Read is not within window => report old and initiate new molecule for that barcode.
                else:
                    all_molecules.report(molecule=molecule)
                    all_molecules.terminate(molecule=molecule)

                    molecule = Molecule(barcode=barcode, start=read_start, stop=read_stop, read_header=read.query_name)
                    all_molecules.cache_dict[molecule.barcode] = molecule

            else:
                molecule = Molecule(barcode=barcode, start=read_start, stop=read_stop, read_header=read.query_name)
                all_molecules.cache_dict[molecule.barcode] = molecule

        elif barcode and read.is_unmapped:
            summary.unmapped_bc_tagged_read += 1

    all_molecules.report_and_remove_all()

    return all_molecules.bc_to_mol, all_molecules.header_to_mol


def fetch_bc(pysam_read, barcode_tag, summary=None):
    """
    Fetches barcode from a bam file tag, returns None if reads isn't tagged.
    """

    try:
        barcode = pysam_read.get_tag(barcode_tag)
    except KeyError:
        barcode = None

        if summary:
            summary.reads_without_barcode += 1

    return barcode


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
        return self.stop - self.start

    def add_read(self, stop, read_header):
        """
        Updates molecule's stop position, number of reads and header name set()
        """

        self.stop = stop
        self.read_headers.add(read_header)

        self.number_of_reads += 1


class AllMolecules:
    """
    Tracks all molecule information, with finished molecules in .bc_to_mol, and molecules which still might get more
    reads in .cache_dict.
    """

    def __init__(self, min_reads):
        """
        :param min_reads: Minimum reads required to add molecule to .bc_to_mol from .cache_dict
        """

        # Min required reads for calling proximal reads a molecule
        self.min_reads = min_reads

        # Molecule tracking system
        self.cache_dict = dict()

        # Dict for finding mols belonging to the same BC
        self.bc_to_mol = dict()

        # Dict for finding mol ID when writing out
        self.header_to_mol = dict()

    def report(self, molecule):
        """
        Commit molecule to .bc_to_mol, if molecule.reads >= min_reads
        """

        if molecule.number_of_reads >= self.min_reads:
            if molecule.barcode not in self.bc_to_mol:
                self.bc_to_mol[molecule.barcode] = set()
            self.bc_to_mol[molecule.barcode].add(molecule)
            for header in molecule.read_headers:
                self.header_to_mol[header] = molecule.ID

    def terminate(self, molecule):
        """
        Removes a specific molecule from .cache_dict
        """

        del self.cache_dict[molecule.barcode]

    def report_and_remove_all(self):
        """
        Commit all .cache_dict molecules to .bc_to_mol and empty .cache_dict (provided they meet criterias by report
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
        self.reads_tagged = int()

        self.overlapping_reads_in_molecule = int()
        self.unmapped_bc_tagged_read = int()
        self.reads_without_barcode = int()
        self.duplicates = int()

    def non_analyzed_reads(self):

        return (
            self.overlapping_reads_in_molecule
            + self.reads_without_barcode
            + self.unmapped_bc_tagged_read
            + self.duplicates)

    def print_stats(self):
        """
        Prints stats to terminal
        """

        logger.info(f"Tot reads in file: {self.tot_reads}")
        logger.info(f"Non-analyzed reads: {self.non_analyzed_reads()}")
        logger.info(f"  Reads overlapping within molecule: {self.overlapping_reads_in_molecule}")
        logger.info(f"  Reads without barcode: {self.reads_without_barcode}")
        logger.info(f"  Unmapped reads with barcode {self.unmapped_bc_tagged_read}")
        logger.info(f"  Duplicate reads: {self.duplicates}")
        logger.info(f"Reads tagged: {self.reads_tagged}")

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
            print(number_of_molecules, file=molecules_per_bc)
            for molecule in (molecule_dict[barcode]):
                print(f"{molecule.number_of_reads}\t{molecule.length()}\t{barcode}\t{number_of_molecules}",
                      file=molecule_stats)

        # Close files
        for output_file in (molecules_per_bc, molecule_stats):
            output_file.close()


def add_arguments(parser):
    parser.add_argument("bam",
                        help="Sorted BAM file tagged with barcode in the same tag as specified in -b/--barcode-tag.")
    parser.add_argument("output",
                        help="Output BAM file with molecule tags found under the tag specified at -m/--molecule-tag "
                             "and molecules number for each barcode under the specified -n/--number-tag.")

    parser.add_argument("-t", "--threshold", type=int, default=4,
                        help="Threshold for how many reads are required for including given molecule in statistics "
                             "(except_reads_per_molecule). Default: %(default)s")
    parser.add_argument("-w", "--window", type=int, default=30000,
                        help="Window size cutoff for maximum distance in between two reads in one molecule. Default: "
                             "%(default)s")
    parser.add_argument("-b", "--barcode-tag", default="BX",
                        help="SAM tag for storing the error corrected barcode. Default: %(default)s")
    parser.add_argument("-s", "--stats-file",
                        help="Write barcode/molecule statistics files.")
    parser.add_argument("-m", "--molecule-tag", default="MI",
                        help="SAM tag for storing molecule index specifying a identified molecule for each barcode. "
                             "Default: %(default)s")
    parser.add_argument("-n", "--number-tag", default="MN",
                        help="SAM tag for storing molecule count for a particular barcode. Default: %(default)s")
