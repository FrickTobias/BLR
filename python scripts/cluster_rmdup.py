#! /usr/bin/env python

def main():


    #
    # Imports & globals
    #

    global args, summaryInstance, output_tagged_bamfile, sys, time, duplicate_position_dict, singleton_duplicate_position
    import pysam, sys, time

    #
    # Argument parsing
    #
    argumentsInstance = readArgs()

    #
    # Initials
    #
    summaryInstance = Summary()

    #
    # Progress
    #
    report_progress('Reading input file and building duplicate position list')

    #
    # Data processing & writing output
    #

    import pdb

    #
    # Find & mark duplicate positions and save paired reads of those positions
    #
    duplicate_position_dict = dict()
    infile = pysam.AlignmentFile(args.input_tagged_bam, 'rb')
    current_read_count = 1000000
    cache_read_tracker = dict()
    cache_readpair_tracker = dict()
    singleton_duplicate_position = dict()
    for read in infile.fetch(until_eof=True):

        # Progress
        summaryInstance.totalReadPairsCount += 1
        if current_read_count < summaryInstance.totalReadPairsCount:
            report_progress("{:,}".format(current_read_count) + ' pairs read')
            current_read_count += 1000000

        # Cache read system
        header = read.query_name
        if header in cache_read_tracker:
            # Fetch mate and remove from tracking dict
            mate = cache_read_tracker[header]
            del cache_read_tracker[header]
        else:
            # Save read for later when mate is found
            cache_read_tracker[header] = read
            continue

        # Save all reads sharing position
        rp_start = mate.get_reference_positions()[0]
        if rp_start in cache_readpair_tracker:
            cache_readpair_tracker[rp_start].append((mate, read))
        else:

            # Send chunk of reads to function
            for the_only_entry in cache_readpair_tracker.values(): process_readpairs(the_only_entry)
            cache_readpair_tracker = dict()
            cache_readpair_tracker[rp_start] = list()
            cache_readpair_tracker[rp_start].append((mate, read))

    sys.stderr.write('YOU MUST CODE DUPLICATE READ SAVING TO DICT AS FUNCTION\n')
    sys.stderr.write('OTHERWISE FIRST READ IN DUPLICATE GROUPS WILL BE MISSED!\n')

    report_progress('Duplicate positions and barcode ID:s from read pairs saved')
    report_progress('Fetching unpaired read duplicate posions & barcode ID:s')

    #
    # For all unpaired reads, save positions for extending duplicate seeds
    #

    cache_position_tracker = dict()
    for unpaired_read in cache_read_tracker.values():

        chromosome = unpaired_read.reference_name
        positions = unpaired_read.get_reference_positions()
        start_stop = (positions[0], positions[-1])


        if not chromosome in cache_position_tracker:
            cache_position_tracker[chromosome] = dict()

        # If start
        if start_stop in cache_position_tracker[chromosome]:
            cache_position_tracker[chromosome][start_stop].append(unpaired_read)

        else:
            duplicate = bool()
            for read in cache_position_tracker[chromosome].values():

                if read.is_duplicate:
                    duplicate = True
                    break

            if duplicate == True:
                for read in cache_position_tracker[chromosome].values():

                    barcode_ID = int(read.get_tag('RG'))

                    # If chr not in dict, add it
                    if not chromosome in singleton_duplicate_position:
                        singleton_duplicate_position[chromosome] = dict()

                    # If this position has no reads yet, add position as key giving empty list as value
                    if not positions in singleton_duplicate_position[chromosome]:
                        singleton_duplicate_position[chromosome][start_stop] = set()

                    # Add read to dictionary
                    singleton_duplicate_position[chromosome][start_stop].add(barcode_ID)

            cache_read_tracker = dict()
            cache_read_tracker[chromosome] = dict()
            cache_read_tracker[chromosome][start_stop] = list()

    infile.close()

    report_progress('Total read count: ' + "{:,}".format(summaryInstance.totalReadPairsCount))

    #
    # Seeding & extending using duplicate readpairs
    #
    duplicates = BarcodeDuplicates()
    report_progress('Seeding barcode ID duplicates')
    for chromosome in duplicate_position_dict:

        for duplicate_start in sorted(duplicate_position_dict[chromosome].keys()):

            # Devide into exactly matching positions (rp1_pos == rp2_pos == rp3_pos...)
            possible_duplicate_seeds = dict()
            for readpair in duplicate_position_dict[chromosome][duplicate_start]:

                mate = readpair[0]
                mate_pos = mate.get_reference_positions()
                mate_start_stop = tuple(mate_pos[::len(mate_pos)-1]) # Get first and last position
                read = readpair[1]
                read_pos = read.get_reference_positions()
                read_start_stop = tuple(read_pos[::len(read_pos)-1]) # Get first and last position
                barcode_ID = int(read.get_tag('RG'))

                if not (mate_start_stop,read_start_stop) in possible_duplicate_seeds:
                    possible_duplicate_seeds[(mate_start_stop, read_start_stop)] = set()
                possible_duplicate_seeds[(mate_start_stop,read_start_stop)].add(barcode_ID)

            # If more than one barcode at specific position is found, seed.
            for possible_seed, barcode_IDs in possible_duplicate_seeds.items():
                if len(barcode_IDs) >= 2:
                    summaryInstance.duplicateSeeds += 1

                    duplicates.seed(barcode_IDs)


    report_progress("{:,}".format(summaryInstance.duplicateSeeds) + ' seeds generated')
    report_progress('Extending seeds using singleton read duplcates')

    #
    # Extending seeds using duplicate singletons
    #
    for chromosome in singleton_duplicate_position:
        for position in singleton_duplicate_position[chromosome]:
            duplicates.extend(list(singleton_duplicate_position[chromosome][position]))

    #
    # Write output
    #

    barcode_ID_merge_dict = duplicates.fetch_significant_seeds()

    # Reduce several step redundancy in merge dict
    report_progress('Reducing several step redundancy in dictionary')
    for barcode_ID_key in sorted(barcode_ID_merge_dict.keys())[::-1]:
        barcode_ID_value = barcode_ID_merge_dict[barcode_ID_key]
        if barcode_ID_value in barcode_ID_merge_dict:
            del barcode_ID_merge_dict[barcode_ID_key]
            barcode_ID_merge_dict[barcode_ID_key] = barcode_ID_merge_dict[barcode_ID_value]

    report_progress('Writing output')
    infile = pysam.AlignmentFile(args.input_tagged_bam, 'rb')
    out = pysam.AlignmentFile(args.output_bam, 'wb', template=infile)

    if args.explicit_merge:
        explicit_merge_file = open(args.explicit_merge, 'w')
        bc_seq_already_written = set()

    for read in infile.fetch(until_eof=True):

        previous_barcode_id = int(read.get_tag('RG'))

        if not previous_barcode_id in barcode_ID_merge_dict:
            out.write(read)
        else:
            new_barcode_id = str(barcode_ID_merge_dict[previous_barcode_id])
            read.set_tag('RG', new_barcode_id, value_type='Z')
            read.query_name = '_'.join(read.query_name.split('_')[:-1]) + '_RG:Z:' + new_barcode_id

            if args.explicit_merge:
                barcode_seq = read.query_name.split()[0].split('_')[-2]
                if barcode_seq in bc_seq_already_written:
                    pass
                else:
                    bc_seq_already_written.add(barcode_seq)
                    explicit_merge_file.write(str(new_barcode_id) + '\t' + str(barcode_seq) + '\t' +str(previous_barcode_id) + '\n')

    report_progress('Analysis finished')

def process_readpairs(list_of_start_stop_tuples):
    """
    Takes readpairs with the same start positions and process them simultaneously (since if one is a duplicate, all
    should be used for seeding duplicates).
    """
    # First look if duplicates are present in the set of readpairs
    mate_dup, read_dup = bool(), bool()
    list_pos = None
    for readpair in list_of_start_stop_tuples:
        mate = readpair[0]
        read = readpair[1]

        # Look if we only find duplicates in read, mate or in whole read pair
        if read.is_duplicate and mate.is_duplicate:
            mate_dup, read_dup = True, True
            break  # Doesn't need to look more, all are True
        elif mate.is_duplicate:
            mate_dup = True
            list_pos = 0
            if read_dup:
                break  # Doesn't need to look more, all are True
        elif read.is_duplicate:
            read_dup = True
            list_pos = 1
            if mate_dup:
                break  # Doesn't need to look more, all are True

    # If duplicates are found, save appropriately dependant if read, mate or both are marked
    for readpair in list_of_start_stop_tuples:
        mate = readpair[0]
        read = readpair[1]
        # Only fetch positions marked as duplicates
        if mate_dup and read_dup:

            # Stats & fetch information
            summaryInstance.totalReadPairsMarkedAsDuplicates += 2
            chromosome = read.reference_name
            start_position = min(read.get_reference_positions())

            # If chr not in dict, add it
            if not chromosome in duplicate_position_dict:
                duplicate_position_dict[chromosome] = dict()

            # Add all chromosomes into singleton dict as not to get KeyError later
            if not chromosome in singleton_duplicate_position:
                singleton_duplicate_position[chromosome] = dict()

            # Fetch chromosome lengths for progress bars later!

            # If this position has no reads yet, add position as key giving empty list as value
            if not start_position in duplicate_position_dict[chromosome]:
                duplicate_position_dict[chromosome][start_position] = []

            # Add mate and read to dictionary for later investigation
            if not start_position in duplicate_position_dict[chromosome]:
                duplicate_position_dict[chromosome][start_position] = list()
            duplicate_position_dict[chromosome][start_position].append((mate, read))

        else:

            # Check if only only one read is a duplicate, then save position & barcode for extending seeds
            if mate_dup or read_dup:

                single_read = (mate, read)[list_pos]

                # Fetch & Stats
                summaryInstance.totalReadPairsMarkedAsDuplicates += 1
                positions = single_read.get_reference_positions()
                chromosome = single_read.reference_name()

                # If chr not in dict, add it
                if not chromosome in singleton_duplicate_position:
                    singleton_duplicate_position[chromosome] = dict()

                # If this position has no reads yet, add position as key giving empty list as value
                if not positions in singleton_duplicate_position:
                    singleton_duplicate_position[chromosome][positions] = list()

                # Add read to dictionary
                singleton_duplicate_position[chromosome][positions].append(int(single_read.get_tag('RG')))


class BarcodeDuplicates(object):
    """
    Tracks barcode ID:s which have readpairs ovelapping with others.
    """

    def __init__(self):
        """
        Initials
        """

        self.seeds = dict()
        self.threshold = args.threshold
        self.seeds_over_threshold = dict()

    def seed(self, barcode_IDs):
        """
        Adds barcode overlap to dictionary with a value of 1. If an overlap is already present, increases value. Give as
        integers.
        """

        # For all barcode IDs
        for barcode_ID in barcode_IDs:
            # Add to dict if not present
            if not barcode_ID in self.seeds: self.seeds[barcode_ID] = dict()
            # Add all other barcodes as values
            for other_barcodes in barcode_IDs:
                # Don't add self or larger barcode ID values
                if other_barcodes >= barcode_ID: continue
                # Increase value with 1 or seed at value 1
                try: self.seeds[barcode_ID][other_barcodes] += 1
                except KeyError:
                    self.seeds[barcode_ID][other_barcodes] = 1

    def extend(self, barcodeIDs):
        """
        Increases value of barcode seeds with 1 for the overlaps provided
        """

        # Only increases value for previously seeded overlaps.
        for barcode_ID in barcodeIDs:
            for other_barcodes in barcode_ID:
                if other_barcodes >= barcode_ID: continue
                if barcode_ID in self.seeds and other_barcodes in self.seeds[barcode_ID]:
                    self.seeds[barcode_ID][other_barcodes] += 1

    def fetch_significant_seeds(self):
        """
        Fetches all barcode overlaps which are above a specified threshold. Default is 0, aka just having been seeded,
        correlating to having shared exact positions between two readpairs with different barcodes.
        """

        # Find the lowest barcode_id value for all barcodes over threshold
        for from_bc_id in self.seeds:
            to_set = set()
            for to_bc_id in self.seeds[from_bc_id]:
                if self.seeds[from_bc_id][to_bc_id] >= self.threshold:
                    to_set.add(to_bc_id)
            # Commit final entry
            if len(to_set) > 0:
                min_bc_id = min(to_set)
                self.seeds_over_threshold[from_bc_id] = min_bc_id

        return self.seeds_over_threshold

def match_clusterid(clusterid_list_one, clusterid_list_two):
    """
    Takes two lists returns matching entries between the two.
    """

    match_set = set()
    for clusterid_one in clusterid_list_one:
        for clusterid_two in clusterid_list_two:

            if clusterid_one == clusterid_two:
                match_set.add(int(clusterid_one))

    match_sorted_list = sorted(match_set)
    return match_sorted_list

def report_matches(current_match_dict, merge_dict, chromosome):#, duplicate_position_list):
    """ Reports matches by updating merge_dict with newfound matches. Dictionary will be mutually exclusive for
    key > value (can't contain key with more than one value)."""


    # Updates merge_dict() with new matches for every position pair and for every match found at those positions.
    for position, duplicate_set in current_match_dict.items():

        duplicate_list = sorted(duplicate_set)
        lower_cluster_id = duplicate_list[0]

        # If multiple matches found on one position, takes one at the time for easier code structure
        for higher_cluster_id in duplicate_list[1:]:
            prev_value = None
            try:
                prev_value = merge_dict[lower_cluster_id]
            except KeyError:
                merge_dict[higher_cluster_id] = lower_cluster_id
            if not prev_value == None:
                if prev_value == lower_cluster_id:
                    pass
                elif prev_value < lower_cluster_id:
                    merge_dict[higher_cluster_id] = prev_value
                    merge_dict[lower_cluster_id] = prev_value
                else:
                    del merge_dict[higher_cluster_id]
                    merge_dict[higher_cluster_id] = lower_cluster_id
                    merge_dict[prev_value] = lower_cluster_id

    return merge_dict

def report_progress(string):
    """
    Writes a time stamp followed by a message (=string) to standard out.
    Input: String
    Output: [date]  string
    """
    sys.stderr.write(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + '\t' + string + '\n')

class ProgressBar(object):
    """
    Writes a progress bar to stderr
    """

    def __init__(self, name, min, max, step):
        # Variables
        self.min = min
        self.max = max
        self.current_position = min
        self.step = step

        # Metadata
        self.two_percent = (self.max-self.min)/50
        self.current_percentage = self.two_percent

        # If two percent, equivalent of one '#', is less than one step length increase the number of # written each step
        if self.two_percent < self.step and not self.max==2:
            self.progress_length = int(50/(self.max-2))
            self.progress_string = '#' * self.progress_length
        elif self.max == 2:
            self.progress_string = '#' * 25
        else:
            self.progress_string = '#'

        # Printing
        report_progress(str(name))
        sys.stderr.write('\n|------------------------------------------------|\n')

    def update(self):
        # If progress is over 2%, write '#' to stdout
        self.current_position += self.step
        if self.current_percentage < self.current_position:
            sys.stderr.write(self.progress_string)
            sys.stderr.flush()
            time.sleep(0.001)
            self.current_percentage += self.two_percent

    def terminate(self):
         sys.stderr.write('\n')

class readArgs(object):
    """ Reads arguments and handles basic error handling like python version control etc."""

    def __init__(self):
        """ Main funcion for overview of what is run. """

        readArgs.parse(self)
        readArgs.pythonVersion(self)

    def parse(self):

        #
        # Imports & globals
        #
        import argparse, multiprocessing
        global args

        parser = argparse.ArgumentParser(description=__doc__)

        # Arguments
        parser.add_argument("input_tagged_bam", help=".bam file tagged with @RG tags and duplicates marked (not taking "
                                                     "cluster id into account).")
        parser.add_argument("output_bam", help=".bam file without cluster duplicates")

        # Options
        parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                           "Not recommended due to different function "
                                                                           "names in python 2 and 3.")
        parser.add_argument("-t", "--threshold", metavar="<INTEGER>", type=int, default=0, help="Threshold for how many additional overlaps "
                                                                            "(other than four exact positions from two "
                                                                            "readpairs) is needed for mergin two barcode "
                                                                            "clusters.")
        parser.add_argument("-e", "--explicit_merge", metavar="<FILENAME>", type=str, help="Writes a file with new_bc_id \\t original_bc_seq")

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
    """ Summarizes chunks"""

    def __init__(self):

        self.totalReadPairsCount = int()
        self.duplicateSeeds = int()
        self.totalReadPairsMarkedAsDuplicates = int()
        self.duplicatePositionWithoutProximity = int()# Duplicates without proximity to other duplicates (=> cannot be cluster duplicate)
        self.readPairsMerged = int() # Rather the count of reads that have changed cluster ID.
        self.ClustersRemovedDueToMerge = int()
        #self.mergeDict = str() # Removed since there are A LOT of clusters being merged
        self.overlap_dict = dict() # Format for tracking # barcode overlap occurrances during writing of file
        self.coupling_dict = dict() # Summarises overlap_dict (bins values)
        self.readable_coupling_dict = str()
        self.log = args.output_bam + '.log'

        with open(self.log, 'w') as openout:
            pass

    def reportMergeDict(self, merge_dict):
        """ Saves a readable string format of merge_dict to write to out."""

        self.ClustersRemovedDueToMerge = len(merge_dict.keys())

    def coupling_analysis(self):
        """ Summarises overlap into a a coupling analysis dictionary (~bins values)."""

        # Iterate over 'droplet defining' clusterids
        for lower_barcode in self.overlap_dict.keys():
            barcodes_in_drop = len(self.overlap_dict[lower_barcode].keys())

            # iterates over translation event and counts # of puts into key according to how many duplicate reads that event had
            for translation_events, number_of_reads in self.overlap_dict[lower_barcode].items():
                try: self.coupling_dict[barcodes_in_drop][number_of_reads] += 1
                except KeyError:
                    try: self.coupling_dict[barcodes_in_drop]
                    except KeyError:
                        self.coupling_dict[barcodes_in_drop] = dict()
                    self.coupling_dict[barcodes_in_drop][number_of_reads] = 1

        self.readable_coupling_dict += 'bc/drop\toverlap\t#occurances\n'
        for number_of_bc_in_droplet in sorted(self.coupling_dict.keys()):

            for number_of_overlapping_reads in sorted(self.coupling_dict[number_of_bc_in_droplet].keys()):

                self.readable_coupling_dict += str(number_of_bc_in_droplet) + '\t' + str(number_of_overlapping_reads) + '\t' + str(self.coupling_dict[number_of_bc_in_droplet][number_of_overlapping_reads]) + '\n'
                self.readable_coupling_dict += str(number_of_bc_in_droplet) + '\t' + str(number_of_overlapping_reads) + '\t' + str(self.coupling_dict[number_of_bc_in_droplet][number_of_overlapping_reads]) + '\n'

    def writeLog(self):

        import time
        with open(self.log, 'a') as openout:
            openout.write(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
            for objectVariable, value in vars(self).items():
                openout.write('\n\n'+str(objectVariable) + '\n' + str(value))

if __name__=="__main__": main()
