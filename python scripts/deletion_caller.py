 #! /usr/bin/env python2

def main():

    #
    # Imports & globals
    #
    global args, summaryInstance, sys, time, pysam, stats, PS_set_H1, PS_set_H2, last_H, outfile, scipy, stats, infile, numpy, qvalue, alpha
    import pysam, sys, time, scipy, numpy, qvalue
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
    bin_results = bins()

    #
    # Analysis start
    #

    # Settings & initals
    num_bins = 1001
    bin_width = args.bin
    alpha = args.alpha

    # Initials: metadata
    window_size = num_bins*bin_width

    # Open bam file
    infile = pysam.AlignmentFile(args.sort_tag_bam, 'rb')

    # Reporting to stderr
    report_progress('Variables set')
    report_progress('Starting indel calling')

    # Loop over chromosomes
    for chromosome_header in infile.header['SQ']:

        chromosome = chromosome_header['SN']
        chromosome_length = chromosome_header['LN']

        # Option for only running on specific chromosomes.
        if args.chromosome:
            if not chromosome == args.chromosome:
                continue

        # Error handling: Check that total length of chromosome is bigger than one window
        if chromosome_length <= window_size:
            report_progress(str(chromosome) + ' is only ' + str(chromosome_length) + ' bp. Cannot call deletions with current \
            bin size. Consider making bin size smaller if this is an important area.')
            report_progress('Skipping ' + str(chromosome))
            continue

        #
        # 1. Initial window
        #

        # Two main variables for tracking values in current window
        bin_list = list()
        bin_pos_list = list(range(0,window_size+bin_width, bin_width))

        # Loop over all bins and count number of unique barcode ID:s found.
        for i in range(len(bin_pos_list) - 1):
            bin_start = bin_pos_list[i]
            bin_stop = bin_pos_list[i + 1]
            num_bc = count_haplotyped_barcodes(chromosome=chromosome, bin_start=bin_start, bin_stop=bin_stop)
            bin_list.append(num_bc)

        # Calculate p value and report to output if significance below threshold (alpha)
        p_value = heterozygous_deletion_test(bin_list)
        summaryInstance.total_number_of_tests[chromosome] = 1

        # Add bin to results
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

            # Fetch the new bins num_bc
            num_bc = count_haplotyped_barcodes(chromosome=chromosome, bin_start=bin_pos_list[-2], bin_stop=bin_pos_list[-1])

            # Shift window for bin value list
            bin_list = bin_list[1:]
            bin_list.append(num_bc)

            # Test for heterozygous deletions
            p_value = heterozygous_deletion_test(bin_list)
            summaryInstance.total_number_of_tests[chromosome] += 1

            # Add bin to results
            bin_results.add_bin(chromosome, bin_pos_list[500], bin_pos_list[501], \
                                num_bc_H1=bin_list[500][1], num_bc_H2=bin_list[500][2], p_value=p_value)
            # Progress bar
            pg_bar.update()

        # Progress bar
        pg_bar.terminate()

    # Close input
    infile.close()

    # Calculate q-values from p-values
    report_progress('Calculating q-values')
    bin_results.adjust_p_values()

    # Write to file (or std out if specified in options)
    report_progress('Writing output')
    if args.STDOUT:
        bin_results.write_to_stdout()
    else:
        with open(args.outfile, 'w') as outfile:
            bin_results.write_as_vcf()

    #
    # Write logfile containing everything in summaryinstance
    #
    summaryInstance.writeToStdErr()

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

    # Two-tailed binomial test to discern whether the middle bin belongs to the normal distribution
    p_value = scipy.stats.binom_test(x=bin_list[500][1:3], p=mean, alternative='two-sided')

    return p_value

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
        p_array = numpy.array(p_values)
        q_array = qvalue.estimate(p_array)

        # build p to q-value dict
        p_to_q_value_dict = dict()
        for i in range(len(q_array)):
            p_to_q_value_dict[p_array[i]] = q_array[i]

        # Adds a q-value at the end of the list for each bin
        for chrom in self.bin_dict.keys():
            for bin in self.bin_dict[chrom].keys():
                pval = self.bin_dict[chrom][bin][2]
                self.bin_dict[chrom][bin].append(p_to_q_value_dict[pval])

    def write_as_vcf(self):
        """
        Writes a vcf format summary to specified output file
        """

        # Header
        outfile.write('CHR\tSTART:STOP\tTYPE\tqval\tpval\t#bcH1:#bcH2\n')
        for chrom in self.bin_dict.keys():
            for positions, result in self.bin_dict[chrom].items():

                # If args.significant - only write significant qvalues to output file
                if args.significant:
                    qval = result[3]
                    if qval <= alpha:
                        outfile.write(str(chrom) + '\t' + str(positions[0]) + ':' + str(positions[1]) + '\tDEL\t')
                        outfile.write(str(result[3]) + '\t' + str(result[2]) + '\t' + str(result[0]) + ':' + str(result[1]) + '\n')
                    else:
                        pass

                # Write everything to output file
                else:
                    outfile.write(str(chrom) + '\t' + str(positions[0]) +  ':' + str(positions[1]) + '\tDEL\t')
                    outfile.write(str(result[3]) + '\t' + str(result[2]) + '\t' + str(result[0]) + ':' + str(result[1]) + '\n')

    def write_to_stdout(self):
        """
        Writes a vcf format summary to stdout.
        """

        # Header
        print('CHR\tSTART:STOP\tTYPE\tqval\tpval\t#bcH1:#bcH2\n'
        for chrom in self.bin_dict.keys():
            for positions, result in self.bin_dict[chrom].items():

                # If args.significant - only write significant qvalues to std out
                if args.significant:
                    qval = result[3]
                    if qval <= alpha:
                        print(str(chrom) + '\t' + str(positions[0]) + ':' + str(positions[1]) + '\tDEL\t' + str(
                            result[3]) + '\t' + str(result[2]) + '\t' + str(result[0]) + ':' + str(result[1]))
                    else:
                        pass

                # Write everything to std out
                else:
                    print(str(chrom) + '\t' + str(positions[0]) + ':' + str(positions[1]) + '\tDEL\t' + str(
                        result[3]) + '\t' + str(result[2]) + '\t' + str(result[0]) + ':' + str(result[1]))

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
        parser.add_argument("-c", "--chromosome", type=str, help="Only run analysis on the specified chromosome.")
        parser.add_argument("-S", "--STDOUT", action="store_true", help="Writes output to stdout instead of outfile. It "
                                                                        "is still needed to specify output file name.")
        parser.add_argument("-s", "--significant", action="store_true", help="Only write significant hits to out. ")

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

    def writeToStdErr(self):
        """
        Writes all object variables to stdout.
        """

        for objectVariable, value in vars(self).items():
            sys.stderr.write('\n\n' + str(objectVariable) + '\n' + str(value))
        sys.stderr.write('\n')

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