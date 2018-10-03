#! /usr/bin python3

def main():


    #
    # Imports & globals
    #

    global args, summaryInstance, output_tagged_bamfile, sys, time, duplicate_position_dict, singleton_duplicate_position, overlapValues, window, duplicates, pos_dict
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
    # Start of script
    #

    report_progress('Starting Analysis')


    #######
    ## 1 ##     Find & mark duplicate positions and save paired reads of those positions
    #######

    duplicate_position_dict = dict()
    infile = pysam.AlignmentFile(args.input_tagged_bam, 'rb')
    cache_read_tracker = dict()
    cache_readpair_tracker = dict()
    singleton_duplicate_position = dict()
    first_read = True
    duplicates = BarcodeDuplicates()
    window = 100000
    pos_dict = dict()
    overlapValues = OverlapValues()
    tot_read_pair_count = int()
    progress = ProgressReporter('Reads processed', 1000000)
    for read in infile.fetch(until_eof=True):

        tot_read_pair_count += 1
        progress.update()

        if first_read:
            prev_chromosome = read.reference_name
            first_read = False

        # Should only use one alignments in calculations, currently assumes primary is correct
        if read.is_secondary:
            summaryInstance.non_primary_alignments += 1
            continue

        # Cache read system
        header = read.query_name
        if header in cache_read_tracker:
            # Fetch mate and remove from tracking dict
            mate = cache_read_tracker[header]
            del cache_read_tracker[header]
            summaryInstance.intact_read_pairs += 1
        else:
            # Save read for later when mate is found
            cache_read_tracker[header] = read
            continue

        # Check if a read or mate is unmapped, if so, send mapped record to unpaired_reads
        if read.is_unmapped and mate.is_unmapped:
            summaryInstance.unmapped_read_pair += 1
            continue
        elif read.is_unmapped:
            summaryInstance.unmapped_singleton += 1
            cache_read_tracker[mate.query_name] = mate
            continue
        elif mate.is_unmapped:
            summaryInstance.unmapped_singleton += 1
            cache_read_tracker[read.query_name] = read
            continue

        # Save all reads sharing position
        all_read_pos = read.get_reference_positions()
        all_mate_pos = mate.get_reference_positions()
        read_start = all_read_pos[0]
        mate_start = all_mate_pos[0]
        read_stop = all_mate_pos[-1]
        mate_stop = all_mate_pos[-1]

        rp_position_tuple = (read_start, read_stop, mate_start, mate_stop)

        if rp_position_tuple in cache_readpair_tracker:
            cache_readpair_tracker[rp_position_tuple].append((mate, read))
        else:

            # Every time a new chromosome is found, send duplicates for processing (seed_duplicates)
            if not read.reference_name == prev_chromosome:
                seed_duplicates(duplicate_position_dict, chromosome=prev_chromosome)
                cache_read_tracker = dict()
                duplicate_position_dict = dict()

            # Send chunk of reads to classification function: two duplicates => duplicate_position_dict
            for the_only_entry in cache_readpair_tracker.values(): process_readpairs(list_of_start_stop_tuples=the_only_entry)
            cache_readpair_tracker = dict()
            cache_readpair_tracker[rp_position_tuple] = list()
            cache_readpair_tracker[rp_position_tuple].append((mate, read))

            prev_chromosome = read.reference_name

    # Takes care of the last chunk of reads
    for the_only_entry in cache_readpair_tracker.values(): process_readpairs(list_of_start_stop_tuples=the_only_entry)
    seed_duplicates(duplicate_position_dict, prev_chromosome)
    duplicate_position_dict = dict()

    report_progress('Total reads in file:\t' + "{:,}".format(tot_read_pair_count))
    report_progress('Total paired reads:\t' + "{:,}".format(summaryInstance.intact_read_pairs*2))
    report_progress('Reads in unmapped read pairs:\t' + "{:,}".format(summaryInstance.unmapped_read_pair*2))
    report_progress('Non-primary alignments in file:\t' + "{:,}".format(summaryInstance.non_primary_alignments))

    # Close input file
    infile.close()

    report_progress('Removing overlaps under threshold and reducing several step redundancy')
    report_progress('Barcodes seeded for removal:\t' + "{:,}".format(len(duplicates.seeds)))

    # Fetch all seeds which are above -t (--threshold, default=0) number of overlaps (require readpair overlap for seed)
    for bc_id_set in duplicates.seeds:
        duplicates.reduce_to_significant_overlaps(bc_id_set)
    report_progress('Barcodes over threshold (' + str(args.threshold) +'):\t' + "{:,}".format(len(duplicates.translation_dict.keys())))

    # Remove several step redundancy (5 -> 3, 3 -> 1) => (5 -> 1, 3 -> 1)
    duplicates.reduce_several_step_redundancy()
    barcode_ID_merge_dict = duplicates.translation_dict
    report_progress('Barcodes removed:\t\t' + "{:,}".format(len(barcode_ID_merge_dict)))
    report_progress('Barcode dict finished')

    # Option: EXPLICIT MERGE - Writes bc_seq + prev_bc_id + new_bc_id
    if args.explicit_merge:
        explicit_merge_file = open(args.explicit_merge, 'w')
        bc_seq_already_written = set()

    # Write output
    progressBar = ProgressBar(name='Writing output', min=0, max=tot_read_pair_count, step=1)
    infile = pysam.AlignmentFile(args.input_tagged_bam, 'rb')
    out = pysam.AlignmentFile(args.output_bam, 'wb', template=infile)
    for read in infile.fetch(until_eof=True):
        previous_barcode_id = int(read.get_tag(args.barcode_tag))

        # If read barcode in merge dict, change tag and header to compensate.
        if not previous_barcode_id in barcode_ID_merge_dict:
            pass
        else:
            new_barcode_id = str(barcode_ID_merge_dict[previous_barcode_id])
            read.set_tag(args.barcode_tag, new_barcode_id, value_type='Z')
            read.query_name = '_'.join(read.query_name.split('_')[:-1]) + '_@' + args.barcode_tag + ':Z:' + new_barcode_id

            # Option: EXPLICIT MERGE - Write bc seq and new + prev bc ID
            if args.explicit_merge:
                barcode_seq = read.query_name.split()[0].split('_')[-2]
                # Only write unique entries.
                if barcode_seq in bc_seq_already_written:
                    pass
                else:
                    bc_seq_already_written.add(barcode_seq)
                    explicit_merge_file.write(str(new_barcode_id) + '\t' + str(barcode_seq) + '\t' +str(previous_barcode_id) + '\n')

        # Write to out and update progress bar
        out.write(read)
        progressBar.update()

    # Terminate progress bar and close output file
    progressBar.terminate()
    out.close()
    if args.explicit_merge: explicit_merge_file.close()
    report_progress('Finished')

def update_cache_dict(pos_dict, chromosome, position, window):
    """
    Shifts position dict to only contain entries withing $window of $position.
    """
    # If just starting over a new chromosome, makes new dict and continues.
    if not chromosome in pos_dict:
        pos_dict[chromosome] = dict()
    else:
        # Sort list from small to big value
        for position_from_dict in sorted(pos_dict[chromosome].keys()):
            # Erase all entries which are not withing $window of $position
            if position_from_dict[-1][-1] + window >= position[0]:
                break
            else:
                del pos_dict[chromosome][position_from_dict]

    return pos_dict

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
                singleton_duplicate_position[chromosome][positions].append(int(single_read.get_tag(args.barcode_tag)))

def process_singleton_reads(chromosome, start_stop, list_of_singleton_reads):
    """
    Saves all reads if one is marked as duplicate since 'original' read will not be marked.
    """

    duplicate = bool()

    # Loop through all reads at current position and look for duplicates
    for read in list_of_singleton_reads:
        if read.is_duplicate:
            duplicate = True
            break

    # If one of the reads were marked as duplicates, save all reads at current position
    if duplicate == True:
        for read in list_of_singleton_reads:

            barcode_ID = int(read.get_tag(args.barcode_tag))

            # If chr not in dict, add it
            if not chromosome in singleton_duplicate_position:
                singleton_duplicate_position[chromosome] = dict()

            # If this position has no reads yet, add position as key giving empty list as value
            if not start_stop in singleton_duplicate_position[chromosome]:
                singleton_duplicate_position[chromosome][start_stop] = set()

            # Add read to dictionary
            singleton_duplicate_position[chromosome][start_stop].add(barcode_ID)

def seed_duplicates(duplicate_position_dict, chromosome):
    """
    Seeds duplicates for read pairs
    :param duplicate_position_dict: Dictionary with read pairs where both read&mate are marked as duplicates
    :return: Nothing, sends results to possible_duplicate_seeds & overlapValues
    """
    global pos_dict

    # If all reads for this chromosome has been "unpaired duplicates", won't be able to seed => return
    if not chromosome in duplicate_position_dict:
        return

    # Function takes long time, so nice to see if something happens
    # UPDATE
    progressBar = ProgressBar(name=str('Seeding ' + str(chromosome)), min=0, max=len(duplicate_position_dict[chromosome]), step=1)

    for duplicate_pos in sorted(duplicate_position_dict[chromosome].keys()):

        # Devide into exactly matching positions (rp1_pos == rp2_pos == rp3_pos...)
        possible_duplicate_seeds = dict()
        for readpair in duplicate_position_dict[chromosome][duplicate_pos]:

            # Fetching information from read/mate
            mate = readpair[0]
            read = readpair[1]
            read_pos_tuple = (read.get_reference_positions()[0], read.get_reference_positions()[-1])
            mate_pos_tuple = (mate.get_reference_positions()[0], mate.get_reference_positions()[-1])
            readpair_pos_tuple = (mate_pos_tuple, read_pos_tuple)
            barcode_ID = int(read.get_tag(args.barcode_tag))

            # Add all barcodes IDs to set, check later if total > 2 at positions.
            if not readpair_pos_tuple in possible_duplicate_seeds:
                possible_duplicate_seeds[readpair_pos_tuple] = set()
            possible_duplicate_seeds[readpair_pos_tuple].add(barcode_ID)

        # When all overlaps found for current position, try finding proximal read pairs with same overlap
        for possible_seed, barcode_IDs in possible_duplicate_seeds.items():  # Only one entry

            # Increase value for duplicate read pair
            overlapValues.add_bc_set(bc_set=barcode_IDs, readpair=True)

            # Update Proximal read dict to only contain read pairs which are close by
            pos_dict = update_cache_dict(pos_dict, chromosome, position=possible_seed[0], window=window)

            for position in pos_dict[chromosome]:

                # Check overlapping/remaining/non-added
                overlapping_bc_ids = barcode_IDs & pos_dict[chromosome][position]

                # Add if more than two are overlapping
                if len(overlapping_bc_ids) >= 2:
                    duplicates.seeds.add(tuple(sorted(overlapping_bc_ids)))

            # Add current set to pos dict for next iteration
            pos_dict[chromosome][readpair_pos_tuple] = barcode_IDs

        # Update for every position parsed
        # UPDATE
        progressBar.update()

    # Terminate progress bar
    # UPDATE
    progressBar.terminate()

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

class OverlapValues(object):
    """
    bc_id^2 matrix tracking overalap values.
    """

    def __init__(self):

        self.matrix = dict()

    def add_bc_set(self, bc_set, readpair):
        """
        Increases overlap value for all bc ids in bc_set
        :param bc_set: Set of barcodes of which to increase overlap values to
        :param readpair: If True, will add value 2 for all overlaps, otherwise adds 1
        :return: -
        """

        for bc_id in bc_set:
            for other_bc_id in bc_set:
                if bc_id < other_bc_id:

                    if not other_bc_id in self.matrix:
                        self.matrix[other_bc_id] = dict()
                        self.matrix[other_bc_id][bc_id] = int()
                    elif not bc_id in self.matrix[other_bc_id]:
                        self.matrix[other_bc_id][bc_id] = int()

                    if readpair:
                        self.matrix[other_bc_id][bc_id] += 2
                    else:
                        self.matrix[other_bc_id][bc_id] += 1

    def fetch_value_vectors(self, bc_set):
        """
        Fetches the maximum overlap value for all barcode ids in the current set
        :param bc_set: Overlap set of seeded barcode ids
        :return: dict with one key/bc_id => overlap value for threshold comparison.
        """

        value_vector_dict = dict()
        for bc_id in bc_set:
            value_vector_dict[bc_id] = list()

            # Fetch all overlap values and save to list for current bc_id.
            for other_bc_id in bc_set:

                if bc_id < other_bc_id:
                    value_vector_dict[bc_id].append(self.matrix[other_bc_id][bc_id])
                elif bc_id > other_bc_id:
                    value_vector_dict[bc_id].append(self.matrix[bc_id][other_bc_id])

            # Remove list and keep max value
            value_vector_dict[bc_id] = max(value_vector_dict[bc_id])

        return value_vector_dict

class BarcodeDuplicates(object):
    """
    Tracks barcode ID:s which have readpairs ovelapping with others.
    """

    def __init__(self):
        """
        Initials
        """

        self.seeds = set()
        self.threshold = args.threshold
        self.translation_dict = dict()

    def reduce_to_significant_overlaps(self, bc_set):
        """
        Fetches all barcode overlaps which are above a specified threshold. Default is 0, aka just having been seeded,
        correlating to having shared exact positions between two readpairs with different barcodes.
        """

        # Fetch max overlap values for the given bc_set
        bc_set = set(bc_set)
        overlap_values = overlapValues.fetch_value_vectors(bc_set)

        # Remove overlaps with value < threshold
        for pot_merge in bc_set.copy():
            if overlap_values[pot_merge] <= self.threshold:
                bc_set.remove(pot_merge)

        # If more than one barcode id has value >= threshold, commit to translation dict
        if len(bc_set) >= 2:
            min_id = min(bc_set)
            for bc_id in bc_set:
                if bc_id > min_id:
                    self.add(bc_id, min_id)

    def reduce_several_step_redundancy(self):
        """
        Takes translation  dict saved in object and makes sure 5->3, 3->1 becomes 5->1, 3->5
        """

        # Goes from high to low value of barcode IDs
        for barcode_to_remove in sorted(self.translation_dict.keys())[::-1]:

            # Fetch value for entry
            barcode_to_keep = self.translation_dict[barcode_to_remove]

            # If value also is key, adjust according to description of function
            if barcode_to_keep in self.translation_dict:
                real_barcode_to_keep = self.translation_dict[barcode_to_keep]
                del self.translation_dict[barcode_to_remove]
                self.translation_dict[barcode_to_remove] = real_barcode_to_keep

    def add(self, bc_id, min_id):
        """
        Adds bc ids to translation dict without introducing multiple values for keys.
        """

        # If key already exist, find out lowest barcode id value
        if bc_id in self.translation_dict:
            other_value = self.translation_dict[bc_id]

            # If it is the same entry, don't do anything
            if other_value == min_id:
                pass
            # If prev entry is lower than min_id, add min_id as key to give prev value
            elif other_value < min_id:
                self.translation_dict[min_id] = other_value
            # If min_id is lower than prev key, adjust value for bc_id and add prev key => min_id
            elif other_value > min_id:
                del self.translation_dict[bc_id]
                self.translation_dict[bc_id] = min_id
                self.translation_dict[other_value] = min_id

        # If not present, just add to dictionary
        else:
            self.translation_dict[bc_id] = min_id

def report_progress(string):
    """
    Writes a time stamp followed by a message (=string) to standard out.
    Input: String
    Output: [date]  string
    """
    sys.stderr.write(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + '\t' + string + '\n')

class ProgressReporter(object):
    """
    Writes to out during iteration of unknown length
    """

    def __init__(self, name_of_process, report_step):

        self.name = name_of_process
        self.report_step = report_step
        self.position = int()
        self.next_limit = report_step

    def update(self):

        self.position += 1
        if self.position >= self.next_limit:
            report_progress(self.name + '\t' + "{:,}".format(self.position))
            self.next_limit += self.report_step

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

class FileReader(object):
    """
    Reads input files, handles gzip.
    """
    def __init__(self, filehandle, filehandle2=None):

        # Init variables setting
        self.filehandle = filehandle
        self.gzip = bool()

        # Open files as zipped or not not (depending on if they end with .gz)
        if self.filehandle[-3:] == '.gz':
            report_progress('File detected as gzipped, unzipping when reading')
            import gzip
            self.openfile = gzip.open(self.filehandle, 'r')
            self.gzip = True
        else:
            self.openfile = open(self.filehandle, 'r')

        # Paired end preparation
        self.filehandle2 = filehandle2
        if self.filehandle2:

            # Open files as zipped or not not (depending on if they end with .gz)
            if self.filehandle2[-3:] == '.gz':
                report_progress('File detected as gzipped, unzipping when reading')
                import gzip
                self.openfile2 = gzip.open(self.filehandle2, 'r')
            else:
                self.openfile2 = open(self.filehandle2, 'r')

    def fileReader(self):
        """
        Reads non-specific files as generator
        :return: lines
        """
        for line in self.openfile:
            if self.gzip:
                line = line.decode("utf-8")
            yield line

    def fastqReader(self):
        """
        Reads lines 4 at the time as generator
        :return: read as fastq object
        """

        line_chunk = list()
        for line in self.openfile:
            if self.gzip:
                line = line.decode("utf-8")
            line_chunk.append(line)
            if len(line_chunk) == 4:
                read = FastqRead(line_chunk)
                line_chunk = list()
                yield read

    def fastqPairedReader(self):
        """
        Reads two paired fastq files and returns a pair of two reads
        :return: read1 read2 as fastq read objects
        """

        line_chunk1 = list()
        line_chunk2 = list()
        for line1, line2 in zip(self.openfile, self.openfile2):
            if self.gzip:
                line1 = line1.decode("utf-8")
                line2 = line2.decode("utf-8")
            line_chunk1.append(line1)
            line_chunk2.append(line2)
            if len(line_chunk1) == 4 and len(line_chunk2) == 4:
                read1 = FastqRead(line_chunk1)
                read2 = FastqRead(line_chunk2)

                # Error handling
                if not read1.header.split()[0] == read2.header.split()[0]:
                    import sys
                    sys.exit('INPUT ERROR: Paired reads headers does not match.\nINPUT ERROR: Read pair number:\t'+str(progress.position+1)+'\nINPUT ERROR: '+str(read1.header)+'\nINPUT ERROR: '+str(read2.header)+'\nINPUT ERROR: Exiting')
                line_chunk1 = list()
                line_chunk2 = list()
                yield read1, read2

    def close(self):
        """
        Closes files properly so they can be re-read if need be.
        :return:
        """
        self.openfile.close()
        if self.filehandle2:
            self.openfile2.close()

class FastqRead(object):
    """
    Stores read as object.
    """

    def __init__(self, fastq_as_line):

        self.header = fastq_as_line[0].strip()
        self.seq = fastq_as_line[1].strip()
        self.comment = fastq_as_line[2].strip()
        self.qual = fastq_as_line[3].strip()

    def fastq_string(self):
        return self.header + '\n' + self.seq  + '\n' + self.comment  + '\n' + self.qual + '\n'

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
        import argparse
        global args

        parser = argparse.ArgumentParser(description="Merges barcodes if they share several read duplicates in proximity "
                                                     "to each other. These 'barcode clusters' arise from several barcodes "
                                                     "being enveloped in the same emulsion.")

        # Arguments
        parser.add_argument("input_tagged_bam", help=".bam file tagged with RG tags and duplicates marked (not taking "
                                                     "cluster id into account).")
        parser.add_argument("output_bam", help=".bam file without cluster duplicates")

        # Options
        parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                           "Not recommended due to different function "
                                                                           "names in python 2 and 3.")
        parser.add_argument("-t", "--threshold", metavar="<INTEGER>", type=int, default=0, help="Threshold for how many additional overlaps "
                                                                            "(other than four exact positions from two "
                                                                            "readpairs) is needed for mergin two barcode "
                                                                            "clusters. DEFAULT: 0")
        parser.add_argument("-e", "--explicit_merge", metavar="<FILENAME>", type=str, help="Writes a file with new_bc_id \\t original_bc_seq")
        parser.add_argument("-bc", "--barcode_tag", metavar="<BARCODE_TAG>", type=str, default='BC', help="Bamfile tag in which the barcode is specified in. DEFAULT: BC")

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
        self.seeded_barcodes = int()
        self.totalReadPairsMarkedAsDuplicates = int()
        self.log = args.output_bam + '.log'
        self.intact_read_pairs = int()
        self.unmapped_read_pair = int()
        self.non_primary_alignments = int()
        self.unmapped_singleton = int()

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
