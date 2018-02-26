#! /usr/bin/env python2

def main():

    #
    # Imports & globals
    #
    global args, summaryInstance, sys, time
    import pysam, sys, time

    #
    # Argument parsing
    #
    argumentsInstance = readArgs()

    #
    # Process data
    #

    report_progress('Fetching a 20 bp barcode from sequence in r1')
    counter = int()
    totalcounter = int()
    limit = 1000000
    # Reads both files at the same time in order to get bc seq into r2 as well
    with open(args.r1) as f1, open(args.r2) as f2:
        with open(args.out_r1, 'w') as openr1, open(args.out_r2, 'w') as openr2:
            for x, y in zip(f1, f2):
                counter += 1
                # Fetch an entire read pair before parsing
                if counter == 1:
                    if not x.split()[0] == y.split()[0]:
                        sys.stderr.write('Headers does not match!')
                        sys.exit()
                    else:
                        header_r1 = x
                        header_r2 = y
                elif counter == 2:
                    seq_r1 = x
                    seq_r2 = y
                elif counter == 3:
                    comment_r1 = x
                    comment_r2 = y
                elif counter == 4:
                    qual_r1 = x
                    qual_r2 = y

                    # Reset
                    counter = 0
                    totalcounter += 1

                    # Readpair fetched, parsing
                    bc_seq = seq_r1[:20]
                    seq_r1_trimmed = seq_r1[20:]
                    qual_r1_trimmed = qual_r1[20:]

                    name_and_pos_r1 = header_r1.split()[0]
                    read_and_index_r1 = header_r1.split()[1]

                    name_and_pos_r2 = header_r2.split()[0]
                    read_and_index_r2 = header_r2.split()[1]

                    openr1.write(name_and_pos_r1 + '_' + bc_seq + ' ' + read_and_index_r1 + '\n' + seq_r1_trimmed + comment_r1 + qual_r1_trimmed)
                    openr2.write(name_and_pos_r2 + '_' + bc_seq + ' ' + read_and_index_r2 + '\n' + seq_r2 + comment_r2 + qual_r2)

                    if totalcounter >= limit:
                        limit += 1000000
                        report_progress(str(totalcounter) + '\treads processed')

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

        # If two percent, equivalent of one '#', is less than one step length increase the number of # written each step
        if self.two_percent < self.step and not self.max==2:
            self.progress_length = int(50/(self.max-2))
            self.progress_string = '#' * self.progress_length
        elif self.max == 2:
            self.progress_string = '#' * 25
        else:
            self.progress_string = '#'

        # Printing
        report_progress(str(name))
        sys.stderr.write('\n|------------------------------------------------|\n')

    def update(self):
        # If progress is over 2%, write '#' to stdout
        self.current_position += self.step
        if self.current_percentage < self.current_position:
            sys.stderr.write(self.progress_string)
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
        parser.add_argument("r1", help="Read 1 fastq file")
        parser.add_argument("r2", help="Read 2 fastq file")
        parser.add_argument("out_r1", help="Read 1 output")
        parser.add_argument("out_r2", help="Read 2 output")

        # Options
        parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                           "Not recommended due to different function "
                                                                           "names in python 2 and 3.")

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

if __name__=="__main__": main()