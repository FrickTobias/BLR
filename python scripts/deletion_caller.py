 #! /usr/bin/env python2

def main():

    #
    # Imports & globals
    #
    global args, summaryInstance, sys, time, pysam, stats, PS_set_H1, PS_set_H2, last_H, outfile, scipy, stats, infile, numpy
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
    #
    # Data processing & writing output
    #

    # Settings & initals
    num_bins = 1001
    bin_width = args.bin
    alpha = args.alpha

    # Initials: metadata
    window_size = num_bins*bin_width

    # Open bam file
    infile = pysam.AlignmentFile(args.sort_tag_bam, 'rb')
    # Open output
    outfile = open(args.outfile, 'w')

    report_progress('Variables set')
    report_progress('Starting indel calling')

    total_p_value_list = list()

    bin_results = bins()

    # Loop over chromosomes
    for chromosome_header in infile.header['SQ']:

        # Error handling: Check that total length of chromosome is bigger than one window
        chromosome = chromosome_header['SN']
        chromosome_length = chromosome_header['LN']
        if chromosome_length <= window_size:
            report_progress(str(chromosome) + ' is only ' + str(chromosome_length) + ' bp. Cannot call deletions with current \
            bin size. Consider making bin size smaller if this is an important area.')
            report_progress('Skipping ' + str(chromosome))
            continue

        # Temporary for only running on chr1
        if not chromosome == 'chr1': break

        #
        # 1. Initial window
        #
        bin_list = list()
        bin_pos_list = list(range(0,window_size+bin_width, bin_width))

        # Loop over all bins and count number of unique barcode ID:s found.
        for i in range(len(bin_pos_list) - 1):
            bin_start = bin_pos_list[i]
            bin_stop = bin_pos_list[i + 1]
            num_bc = count_haplotyped_barcodes(chromosome=chromosome, bin_start=bin_start, bin_stop=bin_stop)
            bin_list.append(num_bc)

        # Calculate p value and report to output if significance below threshold (alpha)
        p_value, mean, raw_value = heterozygous_deletion_test(bin_list)
        total_p_value_list.append(p_value)

        # Summary reporting
        summaryInstance.total_number_of_tests[chromosome] = 1
        if p_value < alpha:
            significantBin = SignificantBin(p_value, mean, bin_pos_list[500], bin_pos_list[501], chromosome, raw_value)
            significantBin.write_outfile()
            summaryInstance.deletions += 1

        bin_results.add_bin(chromosome, bin_pos_list[500], bin_pos_list[501],num_bc_H1=bin_list[500][1], num_bc_H2=bin_list[500][2],p_value=p_value)

        #
        # 2. Main analysis: Counting barcodes and calling heterozygous deletion for the whole chromosome
        #

        # Divide chromosome into bins
        total_bin_pos_list = list(range(num_bins * bin_width, chromosome_length, bin_width))

        # Intiates progress bar
        pg_bar = ProgressBar(name=chromosome, min=total_bin_pos_list[0], max=total_bin_pos_list[-1], step=bin_width)

        # Iterate through all bins
        for new_pos in total_bin_pos_list:

            # Shift window in position list
            bin_pos_list = bin_pos_list[1:]
            bin_pos_list.append(new_pos)

            # Fetch new bin's num_bc
            num_bc = count_haplotyped_barcodes(chromosome=chromosome, bin_start=bin_pos_list[-2], bin_stop=bin_pos_list[-1])

            # Shift window for bin value list
            bin_list = bin_list[1:]
            bin_list.append(num_bc)

            # Test for heterozygous deletions
            p_value, mean, raw_value = heterozygous_deletion_test(bin_list)

            # Track all p values found for q-value conversion
            total_p_value_list.append(p_value)

            # Summary reporting
            summaryInstance.total_number_of_tests[chromosome] += 1

            if p_value < alpha:
                significantBin = SignificantBin(p_value, mean, bin_pos_list[500], bin_pos_list[501], chromosome, raw_value)
                significantBin.write_outfile()
                summaryInstance.deletions += 1

            bin_results.add_bin(chromosome, bin_pos_list[500], bin_pos_list[501], \
                                num_bc_H1=bin_list[500][1], num_bc_H2=bin_list[500][2], p_value=p_value)
            # Progress
            pg_bar.update()
        pg_bar.terminate()

    bin_results.adjust_p_values()
    bin_results.write_to_stdout()

    sys.stdout.write(q_values)


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

        # Summary
        summaryInstance.total_reads += 1

        # Skips read if no barcode tag is present
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

    # Create array for normal distribution numbers
    tmp_bin_list = bin_list[0:475] + bin_list[526:1000]
    distribution_bin_list = [fraction[0] for fraction in tmp_bin_list]

    # Build array for mean and std deviation calculation
    distribution_array = numpy.array(distribution_bin_list)

    # Calculate mean
    mean = numpy.mean(distribution_array)

    # Saving raw data for classification
    raw_value = bin_list[500][0]

    # Two-tailed binomial test to discern whether the middle bin belongs to the normal distribution
    p_value = scipy.stats.binom_test(x=bin_list[500][1:3], p=mean, alternative='two-sided')

    return p_value, mean, raw_value

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

        # Printing
        sys.stderr.write('\n' + str(name))
        sys.stderr.write('\n|------------------------------------------------|\n')

    def update(self):

        # If progress is over 2%, write '#' to stdout
        self.current_position += self.step
        if self.current_percentage < self.current_position:
            sys.stderr.write('#')
            sys.stderr.flush()
            time.sleep(0.001)
            self.current_percentage += self.two_percent

    def terminate(self):

         sys.stderr.write('\n')

class bins(object):
    """
    Tracks all bins with their positions, p-values and q-values
    """
    def __init__(self):

        # Dictionary with key:(start, stop) => value [num_bc_H1, num_bc_H2, p-value, q-value]
        self.bin_dict = dict()

    def add_bin(self, chrom, start, stop, num_bc_H1, num_bc_H2, p_value):

        try: self.bin_dict[chrom]
        except KeyError:
            self.bin_dict[chrom] = dict()
        self.bin_dict[chrom][(start, stop)] = [num_bc_H1, num_bc_H2, p_value]

    def adjust_p_values(self):
        """
        Calculates q-values for all p-values.
        """

        # Fetch all p-values
        p_values = list()
        for chrom in self.bin_dict.keys():
            for bin in self.bin_dict[chrom].keys():
                p_values.append(self.bin_dict[chrom][bin][2])

        # convert to numpy array and calculate q-values
        p_array = numpy.array(self.p_values)
        q_array = qvalues.estimate(p_array)

        # build p to q-value dict
        p_to_q_value_dict = dict()
        for i in range(len(q_values)):
            p_to_q_value_dict[p_array] = q_array[i]

        # Adds a q-value at the end of the list for each bin
        for chrom in self.bin_dict.keys():
            for bin in self.bin_dict[chrom].keys():
                pval = self.bin[chrom][bin][2]
                self.bin[chrom][bin][bin].append(p_to_q_value_dict[pval])

    def write_to_stdout(self):
        """
        Writes a vcf format summary to stdout.
        """
        for chrom in self.bin_dict.keys():
            print(chrom)
            for positions, result in self.bin_dict[chrom].items():
                print(str(positions) + '\t' + str(result))
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

        self.vcf_file_string = str(self.chromosome) + '\t' + str(self.start) + ':' + str(self.stop) + '\tDEL\t' + str(self.p_value) + '\n'
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
        self.total_number_of_tests = dict()
        self.unassigned_reads = int()
        self.missing_bc_seq = int()
        self.total_reads = int()

    def writeToStdOut(self):
        """
        Writes all object variables to stdout.
        """

        for objectVariable, value in vars(self).items():
            sys.stdout.write('\n\n' + str(objectVariable) + '\n' + str(value))
        sys.stdout.write('\n')

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