#! /usr/bin/env python2

def main():


    #
    # Imports & globals
    #
    global args, summaryInstance, output_tagged_bamfile, sys
    import pysam

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

    # Generate dictionaries for every vcf files with deletion positions
    # One haploid and one diploid
    haploid_dict = dict()
    diploid_dict = dict()
    for vcf_file in args.vcf_files:

        haploid_dict[vcf_file] = set()
        diploid_dict[vcf_file] = set()
        vcf_in = pysam.VariantFile(vcf_file)

        for record in vcf_file.fetch(until_eof=True):
            if record.info['DEL']:

                # Disctrete between haploid and diploid
                haploid = True
                diploid = True

                if haploid:
                    # Make dict which gives start:stop (and REMEMBER TO SORT THEM!)
                    haploid_dict[vcf_file].add(record.pos)
                elif diploid = True:
                    diploid_dict[vcf_file].add(record.pos)
                else:
                    # This should be removed after debugging
                    print('NOT HAPLOID NOR DIPLOID! :(\nSomething is wring at with classification.')

        # concatenate overlapping positions so no position is covered twice.
        concatenate_overlapping_positions(start_to_stop_dict=haploid_dict)
        concatenate_overlapping_positions(start_to_stop_dict=diploid_dict)


    infile = pysam.VariantFile(args.input_tagged_bam, 'rb')

    for read in infile.fetch(until_eof=True):
        print(read)
    infile.close()

    # Calculate coverage over for all deletions in bam file

def concatenate_overlapping_positions(start_to_stop_dict):
    """ Merges overlapping start/stop positinos, aka {3:10, 7:15} => {3:15}."""

    for vcf_file in start_to_stop_dict:
        sorted_list_start_pos = sorted(start_to_stop_dict[vcf_file].keys()):

        for i in range(len(sorted_list_start_pos)-1)
            first_start = sorted_list_start_pos[i]
            first_stop = start_to_stop_dict[vcf_file][sorted_list_start_pos[i]]
            j = i +1

            while True:
                second_start = sorted_list_start_pos[j]
                second_stop = start_to_stop_dict[vcf_file][sorted_list_start_pos[j]]

                if first_stop < second_start:
                    # All good
                    break

                else:
                    position_list = sorted([first_start, first_stop, second_start, second_stop])

                    del start_to_stop_dict[vcf_file][first_start]
                    del start_to_stop_dict[vcf_file][second_start]
                    start_to_stop_dict[vcf_file][position_list[0]] = position_list[-1]
                    j += 1

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
        parser.add_argument("input_bam", help=".bam file to check coverage in")
        parser.add_argument("vcf_files", help="Any number of vcf files containing deletion calls", nargs='*')

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

        with open(self.log, 'w') as openout:
            pass

    def reportMergeDict(self, merge_dict):
        """ Saves a readable string format of merge_dict to write to out."""

    def writeLog(self):

        import time
        with open(self.log, 'a') as openout:
            openout.write(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + '\n')
            for objectVariable, value in vars(self).items():
                print('\n'+objectVariable)
                print(value)
                #print(self.attributes)

if __name__=="__main__": main()