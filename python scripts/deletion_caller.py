 #! /usr/bin/env python2

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

    # Settings & initals
    num_bins = 1001
    bin_width = args.bin
    bin_list = list()
    alpha = args.alpha

    # Initials: metadata
    window_size = num_bins*bin_width
    bin_pos_list = list(range(0, window_size, bin_width))

    # Open bam file
    infile = pysam.AlignmentFile(args.sort_tag_bam, 'rb')
    # Open output
    open(args.outfile, 'w')

    # Loop over chromosomes
    for chromosome_header in infile.header['SQ']:

        # Error handling: Check that total length of chromosome is bigger than one window
        chromosome = chromosome_header['SN']
        chromosome_length = chromosome_header['LN']
        if chromosome_length <= window_size:
            report_progress(str(chromosome) + ' is only ' + chromosome_length + ' bp. Cannot call deletions with current \
            bin size. Consider making bin size smaller if this is an important area.')
            report_progress('Skipping ' + str(chromosome))
            continue

        #
        # 1. Initial window
        #

        # Progress
        report_progress('Counting barcodes in initial bins for ' + str(chromosome))

        # Loop over all bins and count number of unique barcode ID:s found.
        for i in range(len(bin_pos_list) - 1):
            bin_start = bin_pos_list[i]
            bin_stop = bin_pos_list[i + 1]
            num_bc = count_haplotyped_barcodes(chromosome=chromosome, bin_start=bin_start, bin_stop=bin_stop)
            bin_list.append(num_bc)

        # Progress
        report_progress('Initial bins counted')
        report_progress('Starting indel calling')

        # Calculate p value and report to output if significance below threshold (alpha)
        p_value, mean, raw_value = heterozygous_deletion_test(bin_list)
        if p_value < alpha:
            significantBin = SignificantBin(p_value, mean, possibility_of_value, bin_start, bin_stop, chromosome, raw_value)
            significantBin.write_outfile()
            summaryInstance.deletions += 1

        #
        # 2. Main analysis: Counting barcodes and calling heterozygous deletion for the whole chromosome
        #

        # Divide chromosome into bins
        total_bin_pos_list = list(range(num_bins * bin_width, chromosome_length, bin_width))

        # Iterate through all bins
        for new_pos in total_bin_pos_list:

            # Shift window in position list
            bin_pos_list = bin_pos_list[1:]
            bin_pos_list.append(new_pos)

            # Fetch new bin's num_bc
            num_bc = count_haplotyped_barcodes(chromosome=chromosome, bin_start=bin_pos_list[-2], bin_stop=bin_pos_list[-1])

            # Shift window for bin value list
            bin_list = bin_list[1:].append(num_bc)

            # Test for heterozygous deletions
            p_value, mean, raw_value = heterozygous_deletion_test(bin_list)
            if p_value < alpha:
                significantBin = SignificantBin(p_value, mean, possibility_of_value, bin_start, bin_stop, raw_value)
                significantBin.write_outfile()
                summaryInstance.deletions += 1

    # Close input
    infile.close()
    # Close output
    outfile.close()

    #
    # Write logfile containing everything in summaryinstance
    #
    summaryInstance.writeToStdOut()

def count_haplotyped_barcodes(chromosome, bin_start, bin_stop):
    """
    Counts how many unique barcodes was found within two positions for each haplotype(bin_pos=tuple(star,stop)).
    Input: Bin coordinates
    Output: tuple(frac_in_H1, num_bc_h1, num_bc_h2)
    """

    # Initials
    barcode_sets = dict()
    barcode_sets[1] = set()
    barcode_sets[2] = set()


    # Fetch all reads within bin
    for read in infile.fetch(str(chromosome), bin_start, bin_stop):

        try: barcode_ID = read.get_tag('BX')
        except KeyError:
            summaryInstance.missing_bc_seq += 1

        # Skips read if haplotype is not assigned.
        try: haplotype = read.get_tag('HP')
        except KeyError:
            summaryInstance.unassigned_reads += 1
            continue

        # Add barcode sequence to set for given barcode
        barcode_sets[haplotype].add(barcode_ID)

    # How many barcodes where assigned to the different Haplotypes
    num_bc_H1 = len(barcode_sets[1])
    num_bc_H2 = len(barcode_sets[2])

    try: H1_H2_distribution = num_bc_H2 / (num_bc_H1 + num_bc_H2)
    except ZeroDivisionError:
        H1_H2_distribution = 0.5
    num_bc = (H1_H2_distribution, num_bc_H1, num_bc_H2)

    return num_bc

def heterozygous_deletion_test(bin_list):
    """
    Calculates p value for h0 = bin500 belongs to the normal distribution fitted from bins[0:475 , 526:1001]
    Input: bin list containing num_bc tuples found in 1001 bins
    Output: p value and corresponding statistical measures
    """
    import numpy

    # Create array for normal distribution numbers
    tmp_bin_list = list()
    print('BIN LIST :::')
    print(len(bin_list))
    print(bin_list[1000])
    print(bin_list[1001])
    print(bin_list[999])
    tmp_bin_list = bin_list[0:475] + bin_list[526:1001]
    distribution_bin_list = [fraction[0] for fraction in tmp_bin_list]

    # Build array for mean and std deviation calculation
    distribution_array = numpy.array(distribution_bin_list)

    # Calculate mean
    mean = numpy.mean(distribution_array)

    # Saving raw data for classification
    raw_value = bin_list[500][0]

    # Two-tailed binomial test to discern whether the middle bin belongs to the normal distribution
    p_value = scipy.stats.binom_test(x=bin_list[500][1:3], p=mean, alternative='two-sided')

    # Summary reporting
    summaryInstance.total_number_of_tests += 1

    return p_value, mean, raw_value

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

    def __init__(self, p_value, mean, bin_start, bin_stop, chromosome, raw_value):

        self.start = bin_start
        self.stop = bin_stop
        self.p_value = p_value
        self.chromosome = chromosome
        self.mean = mean
        self.raw_value = raw_value

    def write_outfile(self):
        """
        Writes all significant heterozygous deletions to output file in vcf format
        """

        self.vcf_file_string = str(self.chromosome) + '\t' + str(self.start) + ':' + str(self.stop) + '\t' + str(self.type) + '\t' + str(self.p_value)
        outfile.write(self.vcf_file_string)

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
        self.unassigned_reads = int()
        self.missing_bc_seq = int()

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