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
        unchecked_positions = list()
        unchecked_positions_set = set()
        pos_to_cluster_id_dict = dict()

        # Continue loop until all duplicates have been saved in merge dict for the current window.
        # Extends window if duplciate is found
        while True:

            new_matches_pos = function(duplicate_position_list, i, j, window)
            pos_to_cluster_id_dict.append(new_matches_pos)

            if len(new_matches_pos) >= 1:
                for match_pos in new_matches_pos:
                    unchecked_positions_set.add(match_pos)
                unchecked_positions = sorted(unchecked_positions_set) # Converting back to list and sorting
                new_matches_pos = []

                ### RETHINK IF REMOVAL WORKS CORRECTLY?
                    # can it be removed and then added again?
                    # Does it matter since merge_dict checks for existing entries?
                    #

                i = unchecked_positions[0]
                unchecked_positions = unchecked_positions[1:]
                unchecked_positions_set.remove(i)
                j = i + 1

            elif len(unchecked_positions) >= 1:
                i = unchecked_positions[0]
                unchecked_positions = sorted(unchecked_positions_set)
                unchecked_positions_set.remove(i)
                j = i + 1

            else:
                true_duplicates = check_duplicate(pos_to_cluster_id_dict, window):    # Take window size into account!
                for duplicate in true_duplicates:

                    lower_cluster_id = duplicate[0]
                    higher_cluster_id = duplicate[1]

                    prev_value = None
                    try: prev_value = merge_dict[lower_cluster_id]
                    except KeyError:
                        merge_dict[higher_cluster_id] = lower_cluster_id
                    if not prev_value == None:
                        if prev_value == lower_cluster_id:
                            pass
                        elif prev_value < cluster_id_low:
                            merge_dict[higher_cluster_id] = prev_value
                            merge_dict[lower_cluster_id] = prev_value
                        else:
                            del merge_dict[higher_cluster_id]
                            merge_dict[higher_cluster_id] = lower_cluster_id
                            merge_dict[prev_value] = lower_cluster_id
                break

    # Merges dict cases where {5:3, 3:1} to {5:1, 3:1} since reads are not ordered according to cluster id.
    for cluster_id_to_merge in sorted(merge_dict.copy.values()):

        # Try to find value for
        try: lower_value = merge_dict[merge_dict[cluster_id_to_merge]]
        except KeyError:
            continue

        higher_value = merge_dict[cluster_id_to_merge]
        del merge_dict[cluster_id_to_merge]
        merge_dict[cluster_id_to_merge] = lower_value
        merge_dict[higher_value] = lower_value

    for read in infile.fetch(until_eof=True):

        # If RG tag i merge dict, change its RG to the lower number
        try: read.tags = str(merge_dict[int(read.tags)])
        except KeyError:
            pass

        # Write read to out.
        out.write(read)

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

def match_bc(barcode_list_one, barcode_list_two):
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