#! /usr/bin/env python2

def main():


    #
    # Imports & globals
    #
    import pysam
    global args, summaryInstance, output_tagged_bamfile

    #
    # Argument parsing
    #
    argumentsInstance = readArgs()
    processor_count = readArgs.processors(argumentsInstance)

    #
    # Initials
    #
    summaryInstance = Summary()

    #
    # Data processing & writing output
    #

    # Consist of three main steps
    #   - Create duplication list/dict
    #   - Filter said list/dict for entries which cannot be duplicates (proximity)
    #   - Check remaining (unique) positions & merge if needed

    ###
    #
    #
    # ANANE'S TIP: pysam.fetch
    # for fetching position specific regions in .bam files (need sort+index)
    #
    #
    ###

    # read duplicate file and build position list which is to be investigated.
    duplicate_position_dict = dict()
    infile = pysam.AlignmentFile(args.input_bam, 'rb')
    for read in infile.fetch(until_eof=True):

        # Only fetch positions marked as duplicates
        if not read.is_duplicate: continue

        # fetch barcode
        cluster_id = read.query_name.split()[0]('_')[-1]
        try: duplicate_position_dict[read.pos].append(cluster_id)
        except KeyError:
            duplicate_position_dict[read.pos] = [cluster_id]

    window = 100000
    proximity_duplication_dict = reduce_dict(duplicate_position_dict, window)  # Window is the strict cutoff value for how far a read can be for it being in the same cluster id.

    duplicate_position_list = sorted(proximity_duplication_dict.copy().keys())
    merge_dict = dict()
    for i in range(len(proximity_duplication_dict.copy().keys()-1)):

        j = i + 1
        while True:

            if duplicate_position_list[i]+window >= duplicate_position_list[j]:
                # if new_matches >= 1
                    # for discovered_match_positions:
                        # unchecked_positions.append(duplicate_position_list[j])
                    # unchecked_positions = sorted(unchecked_positions)
                    # i = unchecked_positions[0]
                    # unchecked_positions = unchecked_positions[1:]
                    # new_matches = 0
                    # j = i + 1
                # elif len(unchecked_position) >= 1 # same as for new_matches >= 1 but skips the adding part.
                    # i = unchecked_positions[0]
                    # unchecked_positions = unchecked_positions[1:]
                    # new_matches = 0
                    # j = i + 1
                # else =>
                    # Check if duplicate criteria met
                        # merge
                        # compare_id = None
                        # try: compare_id = merge_dict[from]
                        # except KeyError:
                            # merge_dict[from] = to
                        # if compare_id:
                            # if compare_id == to:
                                # pass
                            # elif compare_id >= to:
                                # del dict[from] (which gives compare_id currently)
                                # dict[from] = to
                                # dict[compare_id] = to
                            # else:
                                # dict[to].append(compare_id)
                    # break

            # Takes two lists and compared how many matches there is in between the two.
            matches = match_bc(duplicate_position_dict[duplicate_position_list[i]],duplicate_position_dict[duplicate_position_list[j]])

            if len(matches) >= 1:
                match_list.append(matches)
                # For match_pair in matches
                    # Fetch cluster_ids
                    # Sort cluster_ids
                    # store match cluster_id_x => cluster_id_y (from => to)
                    # Number of total matches += 1
                    # New_matches += 1
                    # discovered_match_positions.append(j)

    # ALL READS:
        # if read.tag in merge_dict:
            # merge

    infile.close()

def reduce_dict(unfiltered_position_list, window):
    """ Removes sorted list elements which are not within the window size from the next/previous entry and formats to
    dict instead of list."""

    add_anyway = False
    filtered_position_dict = dict()

    for i in range(len(unfiltered_position_list)-1):
        position=unfiltered_position_list[i]
        next_pos=unfiltered_position_list[i+1]

        if position == next_pos: # The other sequence, fetches sequences for values later
            pass

        elif (position+window) >= next_pos:
            filtered_position_dict.append(position)
            add_anyway = True

        elif add_anyway:
            filtered_position_dict.append(position)
            add_anyway = False

        else:
            pass

    if add_anyway:
        filtered_position_dict.append(position)

    #
    # Stats
    #
    summaryInstance.proximal_duplicates = len(filtered_position_dict.keys())
    summaryInstance.total_positions_marked_as_duplicates = len(unfiltered_position_list)

    return filtered_position_dict

def match_bc(barcode_list_one, barode_list_two):
    """ Takes two lists and compared how many matches there is in between the two. """

    for barcode_one in barcode_list_one:
        for barcode_two in barcode_list_two:
            if barcode_one == barcode_two:
                match_list.append(barcode_one)

    return match_list

class ClusterObject(object):
    """ Cluster object"""

    def __init__(self, clusterId):

        self.barcode_to_bc_dict = dict()
        self.Id = int(clusterId.split()[1]) # Remove 'Cluster' string and \n from end

    def addRead(self, line):

        accession = line.split()[2].rstrip('.')
        barcode = accession.split(':')[-1]
        self.barcode_to_bc_dict[barcode] = self.Id # Extract header and remove '...'

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
        parser.add_argument("-p", "--processors", type=int, default=multiprocessing.cpu_count(),
                            help="Thread analysis in p number of processors. Example: python "
                                 "TagGD_prep.py -p 2 insert_r1.fq unique.fa")

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

    def processors(self):

        #
        # Processors
        #
        import multiprocessing
        processor_count = args.processors
        max_processor_count = multiprocessing.cpu_count()
        if processor_count == max_processor_count:
            pass
        elif processor_count > max_processor_count:
            sys.stderr.write(
                'Computer does not have ' + str(processor_count) + ' processors, running with default (' + str(
                    max_processor_count) + ')\n')
            processor_count = max_processor_count
        else:
            sys.stderr.write('Running with ' + str(processor_count) + ' processors.\n')

        return processor_count

class Summary(object):
    """ Summarizes chunks"""

    def __init__(self):

        #
        # Overall stats
        #
        self.total_positions_marked_as_duplicates = int()
        self.not_real_duplicates = int() # tracks how many duplicate positions which are NOT merged

        #
        # Individual steps
        #
        self.non_proximal_duplicates = int()
        self.proximal_duplicates = int()
        self.proximal_duplicates_with_same_bc_sequence = int()

        self.merged_clusters = int()

        self.log = '.'.join(log) + '.log'
        with open(self.log, 'w') as openout:
            pass

    def updateReadToClusterDict(self, input_dict):
        """ Merges cluster specific dictionaries to a master dictionary."""

        self.CurrentClusterId += 1
        for barcode in input_dict.keys():
            self.read_to_barcode_dict[barcode] = self.CurrentClusterId

    def writeLog(self, line):

        import time
        with open(self.log, 'a') as openout:
            openout.write(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + '\n')
            openout.write(line + '\n')

if __name__=="__main__": main()