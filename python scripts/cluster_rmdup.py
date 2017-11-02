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


    #######
    #
    #   read dup. file and build list
    #
    #   for position in dup_list:
    #
    #       while True:
    #
    #           possible_cluster_dup = find_close_dup()
    #           if possible_cluster_dup==True:
    #
    #               fetch_window_from_bam()
    #               fetch_reads_with_same_bc_in_window()
    #               merge()
    #               ??? if len(uncheck_pos) > 0; don't break
    #               ??? if len(uncheck_pos) i -= 1 (counter on multiplex position stuff)
    #
    #######

    # read duplicate file and build position list which is to be investigated.
    duplicate_position_dict = dict()
    infile = pysam.AlignmentFile(args.input_duplicate_bam, 'rb')
    for read in infile.fetch(until_eof=True):
        try: duplicate_position_list[read.pos] += 1
        except KeyError:
            duplicate_position_dict[read.pos] = 1
    infile.close()

    for duplicate_position in sorted(duplicate_position_dict.copy.keys()):

        position_duplicate = proximal_position(duplicate_position, max_phase_window)

        if position_duplicate:

            read_list = fetch_single_position_reads(position_duplicate[0])
            read_list.append(fetch_single_position_reads(position_duplicate[1]))

            match_bc(read_list)

            # if match barcode:

                positions_updated = True
                while positions_updated:

                    # Fetch ALL reads within phase window
                    for read_to_check in fetched_reads():

                        if match_bc(read_to_check):
                            positions_updated
                            read_to_check.set_tag('@RG', )

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
        parser.add_argument("input_tagged_bam", help=".bam file tagged with @RG tags.")
        parser.add_argument("input_duplicate_bam", help=".bam file with positions marked as duplicates (without regard "
                                                        "to @RG flag.)")
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
        self.read_to_barcode_dict = dict()
        self.CurrentClusterId = 0
        self.barcodeLength = int()
        log = args.output_tagged_bam.split('.')[:-1]
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