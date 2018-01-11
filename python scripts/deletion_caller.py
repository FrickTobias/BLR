#! /usr/bin/env python2

def main():

    #
    # Imports & globals
    #
    global args, summaryInstance, output_tagged_bamfile, sys, time, pysam, stats
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

    # Settings & initals
    num_bins = 1001
    bin_width = args.bin
    bin_list = list()

    # Initials: metadata
    window_size = num_bins*bin_width
    bin_pos_list = range(0, window_size, bin_width)

    # Open bam file
    infile = pysam.AlignmentFile(args.sort_tag_bam, 'rb')

    # Loop over chromosomes
    for chromosome in infile.header['@SQ']:

        # Error handling: Check that total length of chromosome is bigger than one window
        chromosome_length = chromosome.split()[1].split(':')[1]
        if chromosome_length <= window_size:
            report_progress(str(chromosome) + ' is only ' + chromosome_length + ' bp. Cannot call indels with current bin size. Consider making bin size smaller if this is an important area.')
            report_progress('Skipping ' + str(chromosome))
            continue
        else:

        #
        # 1. Initial window
        #

        # Progress
        report_progress('Counting barcodes in initial bins for ' + str(chromosome))

        # Loop over all bins and count number of unique barcode ID:s found.
        for i in range(len(bin_pos_list) - 1):
            bin_start = bin_pos_list[i]
            bin_stop = bin_pos_list[i + 1]
            num_bc = count_haplotyped_barcodes(chromosome=chr1, bin_start=bin_start, bin_stop=bin_stop)
            bin_list.append(num_bc)

        # Progress
        report_progress('Initial bins counted')
        report_progress('Starting indel calling')

        # Test initial window
        # if significant: report

        #
        # 2. Main analysis: calling indels for the whole chromosome
        #

        # Continues until end of chromosome
        # Length of every chromosome given in the third element, but the first is probably only dict part.
        # @SQ   SN:chr1   LN:240000000    ...

        # Divide chromosome into bins
        total_bin_pos_list = range(num_bins * bin_width, chromosome_length, bin_width)

        for new_pos in total_bin_pos_list:

            # Shift window
            bin_pos_list = bin_pos_list[1:].append(new_pos)

            # Fetch new bin's num_bc
            num_bc = count_haplotyped_barcodes(chromosome=chromsome, bin_start=bin_pos_list[-2], bin_stop=bin_pos_list[-1])

            # Shift window for bin value list
            bin_list = bin_list[1:].append(num_bc)

            # Test
            significant = test_bin(bin_list)

            # If significant
            if significant:
                # report significant hit, aka write to out!
                SignificantBin(bin_pos=, p_value=, num_bc=, chromosome=, distribution_mean=)

    # Close input
    infile.close()
    # Close output

    #
    # Write logfile containing everything in summaryinstance
    #
    summaryInstance.writeLog()

def count_haplotyped_barcodes(chromosome, bin_start, bin_stop):
    """
    Counts how many unique barcodes was found within two positions (bin_pos=tuple(star,stop)).
    Input: Bin coordinates
    Output: Count of how many unique barcodes was found within the bin

    """

    # Fetches all reads within window and count unique barcode ID:s found.
    barcode_set = set()

    # Keeps track of which which barcode ID:s which has already been accounted for
    barcode_set_H1 = set()
    barcode_set_H2 = set()

    # Keeps track of which PS ID:s which already have been assigned to a haplotype
    PS_set_H1 = set()
    PS_set_H2 = set()

    for read in infile.fetch(chromosome, bin_start, bin_stop):
        barcode_ID = read.tag['@RG']
        barcode_set.add(barcode_ID)
    num_bc = len(barcode_set)
    return num_bc

def test_bin(bin_list):
    """
    Calculates p-value for h0 = bin500 belongs to the normal distribution fitted from bins[0:475 , 526:1001]
    Input: bin list containing num_bc found in 1001 bins
    Output: p-value and corresponding statistical measures
    """

    # Create array for normal distribution numbers
    normal_distribution_list = bin_list[0:475] + bin_list[526:1001]
    normal_distribution_array = numpy.array(normal_distribution_list)

    # Calculate mean(mu) and standard deviation (sigma)
    std_devitation = numpy.std(normal_distribution_array)
    mean = numpy.mean(normal_distribution_array)

    # Calculates the possibility of to get the number of barcodes that is present in the middle bin, bin500
    z1 = ((bin_list[500]-1)-mean)/std_devitation
    z2 = (bin_list[500]-mean)/std_devitation
    possibility_of_value = z2 - z1

    # Two-tailed binomial test to discern whether the middle bin belongs to the normal distribution
    p-value = scipy.stats.binom_test(x=bin_list[500], n=1, p=possibility_of_value, alternative='two-sided')

    return p-value, normal_distribution_list, mean, std_devitation, possibility_of_value

def report_progress(string):
    """
    Writes a time stamp followed by a message (=string) to standard out.
    Input: String
    Output: [date]  string
    """
    sys.stderr.write(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + '\t' + string + '\n')


class SignificantBin(object):
    """
    Object for saving data for significant deletions.
    """

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
        """
        Writes all significant indels to output file in vcf format
        """

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
        parser.add_argument("output_vcf", help=".vcf file with statistically significantly different barcode abundances.")

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