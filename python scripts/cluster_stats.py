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
    window = args.window_size
    past_unpaired_reads = dict()
    current_phase_block = dict()
    rp_max_dist = 50000

    with pysam.AlignmentFile(args.sort_tag_bam, 'rb') as infile:

        for chromosome_header in infile.header['SQ']:

            chromosome_name = chromosome = chromosome_header['SN']
            chromosome_length = chromosome_header['LN']

            for read in infile.fetch(chromosome_name, 0, chromosome_length):

                summaryInstance.reads += 1

                # Stores read until its mate occurs in file
                # NB: mate will always be upstream of read!
                try: mate = past_unpaired_reads[read.query_name]
                except KeyError:
                    past_unpaired_reads[read.query_name] = read
                    continue

                # Drop read from tmp storage when mate is found
                del past_unpaired_reads[read.query_name]

                mate_start, mate_stop = mate.get_reference_positions()[0], mate.get_reference_positions()[-1]
                read_start, read_stop = read.get_reference_positions()[0], read.get_reference_positions()[-1]

                # Standardise read positions according to refernce instead of read direction.
                mate_start, mate_stop = direct_read_pairs_to_ref(mate_start, mate_stop)
                read_start, read_stop = direct_read_pairs_to_ref(read_stop, read_start)

                # If read pairs are ridiculously far apart, just count reads as
                if mate_stop + rp_max_dist < read_start:
                    summaryInstance.unpaired_reads += 2
                    summaryInstance.unpaired_reads_in_same_chr += 2
                    continue

                percent_coverage = read_pair_coverage(mate_start, mate_stop, read_start, read_stop)

                barcode_id = read.get_tag('RG')

                try: last_pos_of_phase_block = current_phase_block[barcode_id]['stop']
                except KeyError:
                    # Initiate new entry with (start, stop, # reads)
                    # Convert to object or function for initiate
                    current_phase_block[barcode_id] = dict()
                    current_phase_block[barcode_id]['start'] = mate_start
                    current_phase_block[barcode_id]['stop'] = read_stop
                    current_phase_block[barcode_id]['coverage'] = percent_coverage
                    current_phase_block[barcode_id]['number_of_reads'] = 1
                    current_phase_block[barcode_id]['insert_bases'] = mate_stop - mate_start
                    current_phase_block[barcode_id]['bases_btw_inserts'] = 0
                    continue

                # If
                if (last_pos_of_phase_block+window) >= mate_start:
                    current_phase_block[barcode_id]['insert_bases'] = mate_stop - mate_start
                    current_phase_block[barcode_id]['bases_btw_inserts'] = current_phase_block['stop'] - mate_start

                    current_phase_block[barcode_id]['stop'] = read_stop
                    current_phase_block[barcode_id]['coverage'] = (percent_coverage + current_phase_block[barcode_id]['coverage'])
                    current_phase_block[barcode_id]['number_of_reads'] += 1

                else:

                    # Normalises average coverage for number of reads when grand total is known.
                    summaryInstance.reportPhaseBlock(current_phase_block[barcode_id], barcode_id)
                    del current_phase_block[barcode_id]

            # Will rinse/count unpaired reads in between chromosomes
            summaryInstance.unpaired_reads += len(past_unpaired_reads.keys())
            past_unpaired_reads = dict()

            # UNPAIRED READS! Count and report as not proper pairs!

            for barcode_id in current_phase_block.copy():
                summaryInstance.reportPhaseBlock(current_phase_block[barcode_id], barcode_id)
                del current_phase_block[barcode_id]

    summaryInstance.writeResultFiles()

    sys.stderr.write('\nReads found (not read pairs) in bam:\t' + "{:,}".format(summaryInstance.reads) + '\n')
    sys.stderr.write('Unpaired reads (removed from analysis):\t' + "{:,}".format(summaryInstance.unpaired_reads) + '\n')
    sys.stderr.write('In the same chromosome:\t' + "{:,}".format(summaryInstance.unpaired_reads_in_same_chr) + '\n')
    sys.stderr.write('(Defined as being ' + "{:,}".format(rp_max_dist) + ' bp apart)\n')

    # report stats in summary dictionary in understandable way - mean? median? list? plot?

    #summaryInstance.writeToStdErr()

def direct_read_pairs_to_ref(read_start, read_stop):
    """
    Reads can be aligned in forward and backwards direction, this puts all start/stops according to reference positions
    """

    if read_start > read_stop:
        return read_stop, read_start
    else:
        return read_start, read_stop


def read_pair_coverage(mate_start, mate_stop, read_start, read_stop):
    """
    This function calculates the percentage of read bases in the insert.
    """

    if mate_stop >= read_start:
        percent_coverage = 1
    else:
        mate_bp = mate_stop - mate_start
        read_bp = read_stop - read_start
        uncovered_bases = read_start - mate_stop
        percent_coverage = (mate_bp + read_bp) / (uncovered_bases + mate_bp + read_bp)

    return percent_coverage

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
        parser.add_argument("output_prefix", help="prefix for results files (.read_pair_coverage, .coupling_efficiency, "
                                                  ".reads_per_molecule, .molecule_per_barcode and .phase_block_lengths")

        # Options
        parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                           "Not recommended due to different function "
                                                                           "names in python 2 and 3.")
        parser.add_argument("-p", "--processors", type=int, default=multiprocessing.cpu_count(),
                            help="Thread analysis in p number of processors. \nDEFAULT: all available")
        parser.add_argument("-w", "--window_size", type=int, default=100000, help="Window size cutoff for maximum distance "
                                                                                  "in between two reads in one phase block. "
                                                                                  "DEFAULT: 100,000 bp")

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

        self.reads = int()
        self.phase_blocks = int()
        self.bc_clusters = int()

        self.phase_block_lengths = dict()
        self.phase_blocks_per_cluster = dict()
        self.reads_per_phase_block = dict()

        self.ave_coverage_phase_block = int()
        self.ave_bases_read_in_read_pair = int()

        self.phase_block_result_dict = dict()

        self.unpaired_reads = int()
        self.unpaired_reads_in_same_chr = int()

    def reportPhaseBlock(self, phase_block, barcode_id):

        start = phase_block['start']
        stop = phase_block['stop']
        length = stop-start
        num_reads = phase_block['number_of_reads']
        ave_read_pair_coverage = phase_block['coverage'] / num_reads

        # Don't want ZeroDivision error
        if phase_block['bases_btw_inserts'] == 0:
            ave_phase_block_cov = 1
        else:
            ave_phase_block_cov = phase_block['insert_bases']/(phase_block['bases_btw_inserts'] + phase_block['insert_bases'])

        # Tries to append to list of tuples, otherwise creates a tuple list as value for given barcode id
        try: self.phase_block_result_dict[barcode_id]
        except KeyError:
            self.phase_block_result_dict[barcode_id] = list()

        self.phase_block_result_dict[barcode_id].append((start, stop, length, num_reads, ave_read_pair_coverage, ave_phase_block_cov))

    def writeResultFiles(self):

        molecules_per_bc_out = open((args.output_prefix + '.molecules_ber_bc'), 'w')
        coupling_out = open((args.output_prefix + '.coupling_efficiency'), 'w')
        ave_read_pair_coverage_out = open((args.output_prefix + '.ave_read_pair_coverage_in_phase_block'), 'w')
        reads_per_phase_block_out = open((args.output_prefix + '.reads_per_molecule'), 'w')
        phase_block_len_out = open((args.output_prefix + '.phase_block_lengths'), 'w')

        for barcode_id in self.phase_block_result_dict.keys():

            molecules_per_bc_out.write(str(len(self.phase_block_result_dict[barcode_id])) + '\n')

            for phase_block in self.phase_block_result_dict[barcode_id]:

                phase_block_len_out.write(str(phase_block[2]) + '\n')
                reads_per_phase_block_out.write(str(phase_block[3]) + '\n')
                ave_read_pair_coverage_out.write(str(phase_block[4]) + '\n')

                # Not interesting if number of reads found are 1
                if not phase_block[3] == 1:

                    coupling_out.write(str(phase_block[5]/0.5) + '\n')

        for output_file in (molecules_per_bc_out, coupling_out, ave_read_pair_coverage_out, reads_per_phase_block_out, phase_block_len_out):
            output_file.close()


    #def writeToStdErr(self):
    #    """
    #    Writes all object variables to stdout.
    #    """
    #
    #    for objectVariable, value in vars(self).items():
    #        sys.stderr.write('\n\n' + str(objectVariable) + '\n' + str(value))
    #    sys.stderr.write('\n')
    #
    #def writeLog(self):
    #    """
    #    Writes all object variables to a log file (outfile.log)
    #    """
    #
    #    self.log = args.sort_tag_bam + '.log'
    #    import time
    #    with open(self.log, 'w') as openout:
    #        openout.write(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
    #        for objectVariable, value in vars(self).items():
    #            openout.write('\n\n'+str(objectVariable) + '\n' + str(value))

if __name__=="__main__": main()