#! /usr/bin python3

def main():

    #
    # Imports & globals
    #
    global args, summaryInstance, BLR
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
    summaryInstance = Summary()
    window = args.window_size
    currentPhaseBlocks = CurrentPhaseBlocks()

    BLR.report_progress('Running analysis with ' + "{:,}".format(args.window_size) + ' bp window size')
    BLR.report_progress('Fetching reads')
    progress = BLR.ProgressReporter('Reads processed', 1000000)

    # Data processing
    with pysam.AlignmentFile(args.x2_bam, 'rb') as infile:

        # Loop over chromosomes
        for chromosome_header in infile.header['SQ']:

            # Necessary for this way of fetching reads
            chromosome_name = chromosome_header['SN']
            chromosome_length = chromosome_header['LN']

            # For all reads (parsed as single reads and not as read pairs)
            for read in infile.fetch(chromosome_name, 0, chromosome_length):

                progress.update()
                summaryInstance.reads += 1

                # Fetches barcode and skips read removed read pair from stats if not present.
                try: barcode_id = read.get_tag(args.barcode_tag)
                except KeyError:
                    summaryInstance.non_tagged_reads += 1
                    continue

                # Unmapped reads won't affect analysis
                if read.is_unmapped:
                    summaryInstance.unmapped_reads += 1
                    continue

                # Fetch positions
                read_start, read_stop = read.get_reference_positions()[0], read.get_reference_positions()[-1]

                # Standardise read positions according to refernce instead of read direction.
                read_start, read_stop = direct_read_pairs_to_ref(read_stop, read_start)

                # Intitates phase blocks with name = barcode_ID, if not that phase block exist, then fetches last pos(phase block)
                if barcode_id in currentPhaseBlocks.dictionary:
                    last_pos_of_phase_block = currentPhaseBlocks.dictionary[barcode_id]['stop']

                    # If last position of phase block is withing window (100 kb) distance of current read, then add this read to phase block.
                    if (last_pos_of_phase_block+window) >= read_start and last_pos_of_phase_block < read_start:
                        currentPhaseBlocks.addReadPairToPhaseBlock(phase_block=barcode_id, rp_start=read_start, rp_stop=read_stop, query_name=read.query_name)

                    # If reads are overlapping...
                    elif last_pos_of_phase_block >= read_start:

                        # ... add it if it is the mate to a previous read in phaseblock
                        if read.query_name in currentPhaseBlocks.dictionary[barcode_id]['reads']:
                            currentPhaseBlocks.addReadPairToPhaseBlock(phase_block=barcode_id, rp_start=read_start, rp_stop=read_stop, query_name=read.query_name)

                        # ... don't count it if it's just overlapping.
                        else:
                            summaryInstance.overlapping_reads_in_pb += 1

                    # If read is outside range of window from last position, then report the old phase block and initate a new one.
                    else:
                        summaryInstance.reportPhaseBlock(currentPhaseBlocks.dictionary[barcode_id], barcode_id)
                        currentPhaseBlocks.terminatePhaseBlock(phase_block=barcode_id)

                        # Initiate new entry with (start, stop, # reads)
                        currentPhaseBlocks.initiatePhaseBlock(name=barcode_id, start=read_start, stop=read_stop, query_name=read.query_name)

                # Add new phase block if bc_id not in dict yet
                else:
                    currentPhaseBlocks.initiatePhaseBlock(name=barcode_id, start=read_start, stop=read_stop,
                                                          query_name=read.query_name)

            # Report phase blocks when switching to new chromosome, as not to confuse positions
            currentPhaseBlocks.commitAndRemoveAll()

    BLR.report_progress('Phase blocks analysed')

    summaryInstance.writeResultFiles()

    BLR.report_progress('\nReads in bam:\t' + "{:,}".format(progress.position))
    BLR.report_progress('Reads without barcode tag:\t' + "{:,}".format(summaryInstance.non_tagged_reads))
    BLR.report_progress('Overlapping reads within phase_block:\t' + "{:,}".format(summaryInstance.overlapping_reads_in_pb))
    BLR.report_progress('\nMolecules identified:\t' + "{:,}".format(summaryInstance.phase_block_counter))
    BLR.report_progress('Molecules over read threshold (' + str(args.threshold) + '):\t' + "{:,}".format(summaryInstance.phase_blocks_over_threshold))
    if args.filter_bam: BLR.report_progress('Molecules removed:\t' + "{:,}".format(summaryInstance.molecules_over_threshold))
    BLR.report_progress('Drops without more molecules than threshold (' + str(args.threshold) + '):\t' + "{:,}".format(summaryInstance.drops_without_molecules_over_threshold) + '\n')

    with open(args.output_prefix + '.lengths_between_readpairs', 'w') as openin:
        for length, number_of_times in summaryInstance.bp_btw_reads.items():
            for i in range(number_of_times):
                openin.write(str(length) + '\n')

    # Writes output bam file if wanted
    if args.filter_bam:

        progressBar_writeBam = BLR.ProgressBar(name='Writing filtered bam file', min=0, max=summaryInstance.reads, step=1)

        openin = pysam.AlignmentFile(args.x2_bam, 'rb')
        openout = pysam.AlignmentFile(args.filter_bam, 'wb', template=openin)
        for read in openin.fetch(until_eof=True):

            # If no bc_id, just write it to out
            try: barcode_id = read.get_tag('RG')
            except KeyError:
                barcode_id = False

            # If too many molecules in cluster, change tag and header of read
            if barcode_id in summaryInstance.barcode_removal_set:
                tmp_header_list = read.query_name.split('_')
                read.query_name = str(tmp_header_list[0]) + '_' + str(tmp_header_list[1]) + '_RG:Z:NNNNN'
                read.set_tag('RG', 'NNNNN', value_type='Z')
                summaryInstance.reads_with_removed_barcode += 1

            # Writ to out & update progress bar
            openout.write(read)
            progressBar_writeBam.update()

        # End progress bar
        progressBar_writeBam.terminate()

        # Close output
        openout.close()
        openin.close()

    if not summaryInstance.reads == 0:
        BLR.report_progress('Reads with barcodes removed:\t' + "{:,}".format((summaryInstance.reads_with_removed_barcode)) + '\t(' + ("%.2f" % ((summaryInstance.reads_with_removed_barcode/summaryInstance.reads)*100) + ' %)'))
    else:
        BLR.report_progress('No mapped reads found in file.')

def direct_read_pairs_to_ref(read_start, read_stop):
    """
    Reads can be aligned in forward and backwards direction, this puts all start/stops according to reference positions
    """

    # if read backwards, turn it the other way
    if read_start > read_stop:
        return read_stop, read_start
    # Otherwise, return it as is
    else:
        return read_start, read_stop

class CurrentPhaseBlocks(object):
    """
    Tmp storage for phase blocks which might still get more reads assigned to them. Basically a dict with customised functions.
    """

    def __init__(self):
        self.dictionary = dict()

    def initiatePhaseBlock(self, name, start, stop, query_name):

        summaryInstance.phase_block_counter += 1

        self.dictionary[name] = dict()
        self.dictionary[name]['start'] = start
        self.dictionary[name]['stop'] = stop
        self.dictionary[name]['number_of_reads'] = 1
        self.dictionary[name]['bases_btw_inserts'] = 0
        self.dictionary[name]['read_bases'] = stop - start
        self.dictionary[name]['reads'] = set()
        self.dictionary[name]['reads'].add(query_name)

    def addReadPairToPhaseBlock(self, phase_block, rp_start, rp_stop, query_name):

        # Tracks distances between read pairs, won't add value if it is mate to read
        if not query_name in self.dictionary[phase_block]['reads']:
            bp_btw_reads = rp_start - self.dictionary[phase_block]['stop']

            # Tracking the different lengths
            if not bp_btw_reads in summaryInstance.bp_btw_reads:
                summaryInstance.bp_btw_reads[bp_btw_reads] = int()
            summaryInstance.bp_btw_reads[bp_btw_reads] += 1

        # Builds new object
        self.dictionary[phase_block]['read_bases'] += rp_stop - rp_start
        self.dictionary[phase_block]['bases_btw_inserts'] += rp_start - self.dictionary[phase_block]['stop']
        self.dictionary[phase_block]['stop'] = rp_stop
        self.dictionary[phase_block]['number_of_reads'] += 1
        self.dictionary[phase_block]['reads'].add(query_name)

    def terminatePhaseBlock(self, phase_block):

        del self.dictionary[phase_block]

    def commitAndRemoveAll(self):

        for phase_block in self.dictionary.copy().keys():
            summaryInstance.reportPhaseBlock(self.dictionary[phase_block], phase_block)
            del self.dictionary[phase_block]

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
        parser.add_argument("-w", "--window_size", metavar='<INTEGER>', type=int, default=30000, help="Window size cutoff for maximum distance "
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

        # Just defining numbers which will be assigned later
        self.reads = int()
        self.phase_blocks = int()
        self.bc_clusters = int()

        self.phase_block_lengths = dict()
        self.phase_blocks_per_cluster = dict()
        self.reads_per_phase_block = dict()

        self.ave_coverage_phase_block = int()
        self.ave_bases_read_in_read_pair = int()

        self.phase_block_result_dict = dict()

        self.unpaired_reads = int()
        self.unpaired_reads_in_same_chr = int()

        self.phase_blocks_over_threshold = int()
        self.phase_block_counter = int()

        self.non_tagged_reads = int()
        self.drops_without_molecules_over_threshold = int()

        self.overlapping_reads_in_pb = int()

        self.barcode_removal_set = set()
        self.reads_with_removed_barcode = int()

        self.molecules_over_threshold = int()
        self.unmapped_reads = int()

        self.bp_btw_reads = dict()

    def reportPhaseBlock(self, phase_block, barcode_id):

        start = phase_block['start']
        stop = phase_block['stop']
        length = stop-start
        num_reads = phase_block['number_of_reads']
        percent_read_bases = phase_block['read_bases']/(phase_block['bases_btw_inserts'] + phase_block['read_bases'])

        # Tries to append to list of tuples, otherwise creates a tuple list as value for given barcode id
        try: self.phase_block_result_dict[barcode_id]
        except KeyError:
            self.phase_block_result_dict[barcode_id] = list()

        # Save in summary dictionary
        self.phase_block_result_dict[barcode_id].append((start, stop, length, num_reads, percent_read_bases))

    def writeResultFiles(self):

        # Opening all files
        molecules_per_bc_out = open((args.output_prefix + '.molecules_per_bc'), 'w')
        percent_read_bases = open((args.output_prefix + '.percent_read_bases'), 'w')
        reads_per_phase_block_out = open((args.output_prefix + '.reads_per_molecule'), 'w')
        phase_block_len_out = open((args.output_prefix + '.phase_block_lengths'), 'w')
        everything = open((args.output_prefix + '.everything'), 'w')
        #everything.write('read_per_pb\tpb_len\tpb_cov\tave_rp_cov\tmol_per_bc')
        progressBar = BLR.ProgressBar(name='Writing stats files', min=0, max = summaryInstance.phase_block_counter, step = 1)


        # Writing outputs
        for barcode_id in self.phase_block_result_dict.keys():

            molecules_in_cluster = 0
            everything_cache_row = list()
            for phase_block in self.phase_block_result_dict[barcode_id]:

                reads_per_phase_block_out.write(str(phase_block[3]) + '\n')

                # Not interesting if number of reads found are 1
                if phase_block[3] >= args.threshold:
                    molecules_in_cluster += 1
                    percent_read_bases.write(str(phase_block[4]))
                    summaryInstance.phase_blocks_over_threshold += 1
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
                        summaryInstance.molecules_over_threshold += molecules_in_cluster
                        self.barcode_removal_set.add(barcode_id)

                # Writes everything file afterwards in chunks (since it needs molecule per droplet)
                for row in everything_cache_row:
                    everything.write(row + '\t' + str(molecules_in_cluster) + '\n')

        # Close files
        for output_file in (molecules_per_bc_out, percent_read_bases, reads_per_phase_block_out, phase_block_len_out):
            output_file.close()

        progressBar.terminate()

if __name__=="__main__": main()