#! /usr/bin python3

def main():

    #
    # Imports & globals
    #
    global args, summaryInstance, sys, time, pysam, stats, PS_set_H1, PS_set_H2, last_H, outfile, scipy, stats, infile
    import pysam, sys, time, scipy, numpy
    from scipy import stats

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
    # Progress
    #
    report_progress('Building inital binning window')

    #
    # Data processing & writing output
    #

    clusters = ClusterDictObject()

    vcf_results = read_and_process_vcf(args.vcf)

    # Open bam file
    infile = pysam.AlignmentFile(args.sort_tag_bam, 'rb')

    for read in infile.fetch(until_eof=True):

        #
        read_haplotype = read_assigner(read, vcf_results)

        clusters.add_read(read, read_haplotype)


    clusters.adjust_conflicts()
    clusters.write_consensus()
    clusters.write_conflicting()


    # Close input
    infile.close()

    #
    # Write logfile containing everything in summaryinstance
    #
    summaryInstance.writeToStdOut()

def read_and_process_vcf():
    """
    Reads a vcf file and saves the different haplotypes defined by SNV:s
    """



class clusterDictObject():
    """
    Tracks all clusters and is able to give reads from them.
    """

    def ___init__(self):

        self.clusters = dict()
        self.

    def add_read(self, read, haplotype):

        # Fetch barcode ID

        # Add read with barcode ID as key

        # Add value to clusters haplotype

    def adjust_conflitcting_haplotypes(self):

        for cluster in self.clusters.keys():

            # Measure how much many reads contradict each other

            # Measure

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
        import argparse, multiprocessing
        global args

        parser = argparse.ArgumentParser(description=__doc__)

        # Arguments
        parser.add_argument("sort_tag_bam", help=".bam file tagged with @RG tags and duplicates marked (not taking "
                                                     "cluster id into account).")
        parser.add_argument("outfile", help=".vcf file with statistically significantly different barcode abundances.")

        # Options
        parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                           "Not recommended due to different function "
                                                                           "names in python 2 and 3.")
        parser.add_argument("-p", "--processors", type=int, default=multiprocessing.cpu_count(),
                            help="Thread analysis in p number of processors. \nDEFAULT: all available")
        parser.add_argument("-b", "--bin", type=int, default=1000, help="Bin width in bp \nDEFAULT: 1000")
        parser.add_argument("-a", "--alpha", type=float, default=0.05, help="Cutoff for what is considered a significant"
                                                                            "p value. \nDEFAULT: 0.05")

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

    def __init__(self):

        self.deletions = int()
        self.total_number_of_tests = int()

    def writeToStdOut(self):
        """
        Writes all object variables to stdout.
        """

        for objectVariable, value in vars(self).items():
            sys.stdout.write('\n\n' + str(objectVariable) + '\n' + str(value))

    def writeLog(self):
        """
        Writes all object variables to a log file (outfile.log)
        """

        self.log = args.outfile + '.log'
        import time
        with open(self.log, 'w') as openout:
            openout.write(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
            for objectVariable, value in vars(self).items():
                openout.write('\n\n'+str(objectVariable) + '\n' + str(value))

if __name__=="__main__": main()