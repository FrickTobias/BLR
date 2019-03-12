#! /usr/bin python3

def main():

    #
    # Imports & globals
    #
    global args, summary, BLR
    import BLR_functions as BLR, sys, pysam

    #
    # Argument parsing
    #
    argumentsInstance = readArgs()

    # Check python3 is being run
    if not BLR.pythonVersion(args.force_run): sys.exit()

    #
    # Initials
    #
    summary = Summary()
    phaseBlocks = PhaseBlocks()
    prev_chrom = 'chr1'

    # Open file, loop over all reads
    BLR.report_progress('Running analysis with ' + "{:,}".format(args.window) + ' bp window size')
    BLR.report_progress('Fetching reads')
    progress = BLR.ProgressReporter('Reads processed', 1000000)
    with pysam.AlignmentFile(args.x2_bam, 'rb') as infile:
        for read in infile.fetch(until_eof=True):

            # Commit phaseBlocks between chromosomes
            if not prev_chrom == read.reference_name:
                phaseBlocks.reportAndRemoveAll()
                prev_chrom = read.reference_name

            # Fetches barcode and genomic position. Position will be formatted so start < stop.
            BC_id, read_start, read_stop = fetch_and_format(read)
            # If read is unmapped or does not have barcode, skip
            if BC_id == None or read_start == 'unmapped': continue

            # If BC_id already has seen prior reads
            if BC_id in phaseBlocks.dictionary:
                phase_block_stop = phaseBlocks.dictionary[BC_id]['stop']

                # Read is within window => add read to phase block.
                if (phase_block_stop+args.window) >= read_start and phase_block_stop < read_start:
                    phaseBlocks.addRead(name=BC_id, read_start=read_start, read_stop=read_stop, read_name=read.query_name)

                # Overlapping reads => If not overlapping to it's mate, discard read.
                elif phase_block_stop >= read_start:
                    if read.query_name in phaseBlocks.dictionary[BC_id]['reads']:
                        phaseBlocks.addRead(name=BC_id, read_start=read_start, read_stop=read_stop, read_name=read.query_name)
                    else:
                        summary.overlapping_reads_in_pb += 1

                # Read is not within window => report old and initate new phase block.
                else:
                    summary.reportPhaseBlock(name=BC_id, phase_block=phaseBlocks.dictionary[BC_id])
                    phaseBlocks.terminate(name=BC_id)
                    phaseBlocks.initiate(name=BC_id, start=read_start, stop=read_stop, read_name=read.query_name)

            # No previous reads for this bc has been discovered
            else:
                phaseBlocks.initiate(name=BC_id, start=read_start, stop=read_stop, read_name=read.query_name)

            # Progress reporting
            progress.update()
            summary.reads += 1

    # Commit last chr phaseBlocks and log stats
    phaseBlocks.reportAndRemoveAll()
    summary.reads = progress.position
    BLR.report_progress('Phase blocks analysed')

    # Stats to output files and stdout
    summary.writeResultFiles()
    summary.printStats()

    with open(args.output_prefix + '.lengths_between_readpairs', 'w') as openout:
        for length, number_of_times in summary.bp_btw_reads.items():
            for i in range(number_of_times):
                openout.write(str(length) + '\n')

    # Writes output bam file if wanted
    if args.filter_bam:

        progressBar_writeBam = BLR.ProgressBar(name='Writing filtered bam file', min=0, max=summary.reads, step=1)

        openin = pysam.AlignmentFile(args.x2_bam, 'rb')
        openout = pysam.AlignmentFile(args.filter_bam, 'wb', template=openin)
        for read in openin.fetch(until_eof=True):

            # If no bc_id, just write it to out
            try: BC_id = read.get_tag('RG')
            except KeyError:
                barcode_id = False

            # If too many molecules in cluster, change tag and header of read
            if barcode_id in summary.barcode_removal_set:
                tmp_header_list = read.query_name.split('_')
                read.query_name = str(tmp_header_list[0]) + '_' + str(tmp_header_list[1]) + '_RG:Z:NNNNN'
                read.set_tag('RG', 'NNNNN', value_type='Z')
                summary.reads_with_removed_barcode += 1

            # Writ to out & update progress bar
            openout.write(read)
            progressBar_writeBam.update()

        # End progress bar
        progressBar_writeBam.terminate()

        # Close output
        openout.close()
        openin.close()

    try: BLR.report_progress('Reads with barcodes removed:\t' + "{:,}".format((summary.reads_with_removed_barcode)) + '\t(' + ("%.2f" % ((summary.reads_with_removed_barcode/summary.reads)*100) + ' %)'))
    except ZeroDivisionError: BLR.report_progress('No mapped reads found in file.')

def fetch_and_format(read):
    """

    :param read:
    :return:
    """

    try: BC_id = read.get_tag(args.barcode_tag)
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

    return BC_id, read_start, read_stop

class PhaseBlocks(object):
    """
    Tmp storage for phase blocks which might still get more reads assigned to them. Basically a dict with customised functions.
    """

    def __init__(self):
        self.dictionary = dict()

    def initiate(self, name, start, stop, read_name):

        summary.molecules += 1

        self.dictionary[name] = dict()
        self.dictionary[name]['start'] = start
        self.dictionary[name]['stop'] = stop
        self.dictionary[name]['number_of_reads'] = 1
        self.dictionary[name]['bases_btw_inserts'] = 0
        self.dictionary[name]['bases_read'] = stop - start
        self.dictionary[name]['reads'] = set()
        self.dictionary[name]['reads'].add(read_name)

    def addRead(self, name, read_start, read_stop, read_name):

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

    def terminate(self, name):

        del self.dictionary[name]

    def reportAndRemoveAll(self):

        for BC_id in self.dictionary.copy().keys():
            summary.reportPhaseBlock(name=BC_id, phase_block=self.dictionary[BC_id])
            del self.dictionary[BC_id]

class readArgs(object):
    """
    Reads arguments and handles basic error handling like python version control etc.
    """

    def __init__(self):

        readArgs.parse(self)
        readArgs.pythonVersion(self)

    def parse(self):

        #
        # Imports & globals
        #
        import argparse
        global args

        parser = argparse.ArgumentParser(description="Removes barcode tags present at more than -M loci (corresponding "
                                                     "to removing barcode tags from reads origin to droplets which had "
                                                     "more than -M molecules).")

        # Arguments
        parser.add_argument("x2_bam", help=".bam file tagged with @RG tags and duplicates removed. Needs to be indexed.")
        parser.add_argument("output_prefix", help="prefix for results files (.coupling_efficiency, "
                                                  ".reads_per_molecule, .molecule_per_barcode and .phase_block_lengths). "
                                                  "Will also write a .everything file containing rows (equivalent to "
                                                  "molecules) with all the information from the other files.")

        # Options
        parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                           "Not recommended due to different function "
                                                                           "names in python 2 and 3. DEFAULT: False")
        parser.add_argument("-t","--threshold", metavar='<INTEGER>', type=int, default=4, help="Threshold for how many reads are required for including given phase block in statistics (except_reads_per_molecule). DEFAULT: 4")
        parser.add_argument("-w", "--window", metavar='<INTEGER>', type=int, default=30000, help="Window size cutoff for maximum distance "
                                                                                  "in between two reads in one phase block. "
                                                                                  "DEFAULT: 30000")
        parser.add_argument("-bc", "--barcode_tag", metavar='<STRING>', type=str, default='BC', help="Bam file tag where barcode is stored. DEFAULT: BC")
        parser.add_argument("-f", "--filter_bam", metavar='<OUTPUT>', type=str, default=False, help="Write an output bam file where reads have their barcode tag removed if it is "
                                                                                "from a cluster with more molecules than -M (--Max_molecules) molecules. DEFAULT: False")
        parser.add_argument("-M", "--Max_molecules", metavar='<INTEGER>', type=int, default=500, help="When using -f (--filter) this will remove barcode tags for those clusters which have more than -M molecules. DEFAULT: 500")
        args = parser.parse_args()

    def pythonVersion(self):
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

class Summary(object):

    def __init__(self):

        # Stats
        self.reads = int()
        self.phase_blocks = int()
        self.phase_blocks_over_threshold = int()
        self.phase_block_result_dict = dict()
        self.molecules = int()
        self.non_tagged_reads = int()
        self.drops_without_molecules_over_threshold = int()
        self.overlapping_reads_in_pb = int()
        self.barcode_removal_set = set()
        self.reads_with_removed_barcode = int()
        self.unmapped_reads = int()
        self.bp_btw_reads = dict()

        # Filter bam file
        self.molecules_in_outbam = int()
        self.bc_in_outbam = int()

    def reportPhaseBlock(self, name, phase_block):

        start = phase_block['start']
        stop = phase_block['stop']
        length = stop-start
        num_reads = phase_block['number_of_reads']
        percent_bases_read = phase_block['bases_read']/(phase_block['bases_btw_inserts'] + phase_block['bases_read'])

        # Tries to append to list of tuples, otherwise creates a tuple list as value for given barcode id
        try: self.phase_block_result_dict[name]
        except KeyError:
            self.phase_block_result_dict[name] = list()

        # Save in summary dictionary
        self.phase_block_result_dict[name].append((start, stop, length, num_reads, percent_bases_read))

    def printStats(self):

        # Read stats
        BLR.report_progress('\nReads total:\t' + "{:,}".format(self.reads))
        BLR.report_progress('Unmapped reads:\t' + "{:,}".format(self.unmapped_reads))
        BLR.report_progress('Reads without ' + args.barcode_tag + ' tag:\t' + "{:,}".format(self.non_tagged_reads))
        BLR.report_progress('Reads overlapping within phase_block:\t' + "{:,}".format(self.overlapping_reads_in_pb))

        # Molecule stats
        BLR.report_progress('\nMolecules total:\t' + "{:,}".format(self.molecules))
        BLR.report_progress('Molecules kept for stats (min read: ' + str(args.threshold) + '):\t' + "{:,}".format(
            self.phase_blocks_over_threshold))
        BLR.report_progress(
            'BC consequently removed:\t' + "{:,}".format(self.drops_without_molecules_over_threshold) + '\n')

        # Filtering stats
        if args.filter_bam:
            BLR.report_progress('\nMolecules in output bam:\t' + "{:,}".format(self.molecules_in_outbam))
            BLR.report_progress('BC in output bam:\t' + "{:,}".format(self.bc_in_outbam))

    def writeResultFiles(self):

        # Opening all files
        molecules_per_bc_out = open((args.output_prefix + '.molecules_per_bc'), 'w')
        percent_bases_read = open((args.output_prefix + '.percent_bases_read'), 'w')
        reads_per_phase_block_out = open((args.output_prefix + '.reads_per_molecule'), 'w')
        phase_block_len_out = open((args.output_prefix + '.phase_block_lengths'), 'w')
        everything = open((args.output_prefix + '.everything'), 'w')
        progressBar = BLR.ProgressBar(name='Writing stats files', min=0, max = summary.molecules, step = 1)

        # Writing outputs
        for barcode_id in self.phase_block_result_dict.keys():

            molecules_in_cluster = 0
            everything_cache_row = list()
            for phase_block in self.phase_block_result_dict[barcode_id]:

                reads_per_phase_block_out.write(str(phase_block[3]) + '\n')

                # Filter molecules with less than N reads (--threshold)
                if phase_block[3] >= args.threshold:
                    molecules_in_cluster += 1
                    percent_bases_read.write(str(phase_block[4]))
                    summary.phase_blocks_over_threshold += 1
                    phase_block_len_out.write(str(phase_block[2]) + '\n')
                    everything_cache_row.append((str(phase_block[3]) + '\t' + str(phase_block[2]) + '\t' + str(phase_block[4]) + '\t' + str(phase_block[4]) + '\t'  + str(barcode_id) + '\t'))

                progressBar.update()

            # Skips clusters which as a result of -t does not have any molecules left.
            if molecules_in_cluster == 0:
                self.drops_without_molecules_over_threshold += 1
            else:
                molecules_per_bc_out.write(str(molecules_in_cluster) + '\t')

                if args.filter_bam:
                    if molecules_in_cluster > args.Max_molecules:
                        summary.bc_in_outbam += 1
                        summary.molecules_in_outbam += molecules_in_cluster
                        self.barcode_removal_set.add(barcode_id)

                # Writes everything file afterwards in chunks (since it needs molecule per droplet)
                for row in everything_cache_row:
                    everything.write(row + '\t' + str(molecules_in_cluster) + '\n')

        # Close files
        for output_file in (molecules_per_bc_out, percent_bases_read, reads_per_phase_block_out, phase_block_len_out):
            output_file.close()

        progressBar.terminate()

if __name__=="__main__": main()