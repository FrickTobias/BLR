#! /usr/bin/env python2

def main():

    #
    # Imports & globals
    #
    import pysam, sys, time
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
    # Progress
    #
    report_progress('Reading input file and building duplicate position list')

    #
    # Data processing & writing output
    #

    # Settings & initals
    num_bins = 1001
    bin_width = args.bin
    bin_list = list()

    # Initials: metadata
    bin_pos_list = range(0, num_bins * bin_width, bin_width)

    # Open bam file
    infile = pysam.AlignmentFile(args.sort_tag_bam, 'rb')

    for i in range(len(bin_temp_list)-1):

        bin_start = bin_pos_list[i]
        bin_stop = bin_pos_list[i+1]
        barcode_set = set()

        for read in infile.fetch(chr1, bin_start, bin_stop):

            barcode_ID = read.tag['@RG']
            barcode_set.add(barcode_ID)

        bin_list.append(len(barcode_set))

    # Test first window

    for chromsome in infile.header['@SQ']:

        # Continues until end of chromosome
        while True:

            # shift window
            window[]

            # Fetch new bin's num_bc

            # Test
            significant = test_bin(bin_list)

            # If significant
            if significant:
                SignificantBin(bin_pos=, p_value=, num_bc=, chromosome=, distribution_mean=)


    # Close input
    infile.close()
    # Close output

    #
    # Write logfile containing everything in summaryinstance
    #
    summaryInstance.writeLog()

def count_barcodes(bin_pos):
    """ Counts how many unique barcodes was found within two positions (bin_pos=tuple(star,stop))"""

    num_bc = int()

    # add bc seq to set()

    # num_bc = len(list(set))

    return num_bc

def test_bin(bin_list):

    significant = bool()
    p-value = 1

    # Create normal distribution from 475 first and last values
    # Calculate distribution_mean

    # Calculate P(num_bc(bin_500), current_normal_distribution)
    #a = np.array([0.7972, 0.0767, 0.4383, 0.7866, 0.8091, 0.1954, 0.6307, 0.6599, 0.1065, 0.0508])
    #from scipy import stats
    #stats.zscore(a)

    # Test bin 500 for significant difference: two-tailed binomial test
    # x = sum(51 middle bins)
    #scipy.stats.binom_test(x, n=51, p=0.5, alternative='two-sided')

    return significant, p-value, current_normal_distribution, distribution_mean

def report_progress(string):
    """ Writes a time stamp followed by a message (=string) to standard out."""
    sys.stderr.write(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + '\t' + string + '\n')


class SignificantBin(object):
    """ Object for saving data for significant deletions."""

    def __init__(self, bin_pos, p_value, num_bc, chromosome, distribution_mean):

        self.bin_pos = bin_pos
        self.p_value = p_value
        self.num_bc = num_bc
        self.chromosome = chromosome

        if self.num_bc > distribution_mean:
            self.duplication = True
            self.deletion = False
        elif self.num_bc < distribution_mean:
            self.duplication = False
            self.deletion = True

    def write_to_out(self):
        """ Writes all significant indels to output file in vcf format"""

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
        """ Main funcion"""

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
        parser.add_argument("vcf", help=".bam file without cluster duplicates")

        # Options
        parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                           "Not recommended due to different function "
                                                                           "names in python 2 and 3.")
        parser.add_argument("-p", "--processors", type=int, default=multiprocessing.cpu_count(),
                            help="Thread analysis in p number of processors. Example: python "
                                 "TagGD_prep.py -p 2 insert_r1.fq unique.fa")
        parser.add_argument("-b", "--bin", type=int, default=1000, help="Bin width (bp)")

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
    """ Summary"""

    def __init__(self):

        self.log = args.output_bam + '.log'

        with open(self.log, 'w') as openout:
            pass

    def writeLog(self):

        import time
        with open(self.log, 'a') as openout:
            openout.write(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
            for objectVariable, value in vars(self).items():
                openout.write('\n\n'+str(objectVariable) + '\n' + str(value))

if __name__=="__main__": main()