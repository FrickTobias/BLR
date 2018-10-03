#! /usr/bin python3

def main():
    """Takes a fastq file barcode sequences in the header and writes a barcode fasta file with only unique entries. """

    #
    # Imports & globals
    #
    global args, summaryInstance, output_tagged_bamfile, sys, time
    import pysam, sys, time, os

    #
    # Argument parsing
    #
    argumentsInstance = readArgs()

    #
    # Data processing & writing output
    #

    # Generate dict with bc => bc_cluster consensus sequence
    report_progress("Starting analysis")
    clstr_generator = FileReader(args.input_clstr)
    cluster_dict = readAndProcessClusters(clstr_generator.fileReader())
    clstr_generator.close()

    # Read bam files and translate bc seq to BC cluster ID + write to out
    progress = ProgressReporter('Reads processed', 1000000)
    infile = pysam.AlignmentFile(args.input_mapped_bam, 'rb')
    out = pysam.AlignmentFile(args.output_tagged_bam, 'wb', template=infile)
    reads_with_non_clustered_bc = int()
    for read in infile.fetch(until_eof=True):
        read_bc = read.query_name.split()[0].split('_')[-1]

        # Fetch barcode cluster ID based on barcode sequence
        if not read_bc in cluster_dict:
            reads_with_non_clustered_bc += 1
        else:
            bc_id = cluster_dict[read_bc]
            read.set_tag('BC', str(bc_id), value_type='Z')  # Stores as string, makes duplicate removal possible. Can do it as integer as well.
            read.query_name = (read.query_name + '_BC:Z:' + str(bc_id))

        out.write(read)
        progress.update()

    infile.close()
    out.close()
    report_progress('Finished')

def readAndProcessClusters(openInfile):
    """
    Reads clstr file and builds bc => bc_cluster dict (bc_cluster is given as the consensus sequence).
    """

    # For first loop
    if args.skip_nonclust: seqs_in_cluster = 2

    # Reads cluster file and saves as dict
    cluster_dict = dict()
    cluster_ID = int()
    for line in openInfile:

        # Reports cluster to master dict and start new cluster instance
        if line.startswith('>'):

            # If non-clustered sequences are to be omitted, removes if only one sequence makes out the cluster
            if args.skip_nonclust and seqs_in_cluster < 2:
                del cluster_dict[current_key]
            seqs_in_cluster = 0

            cluster_ID += 1
            current_value = cluster_ID
        else:
            current_key = line.split()[2].lstrip('>').rstrip('...').split(':')[2]
            cluster_dict[current_key] = current_value
            seqs_in_cluster +=1

    return(cluster_dict)

def report_progress(string):
    """
    Writes a time stamp followed by a message (=string) to standard out.
    Input: String
    Output: [date]  string
    """
    sys.stderr.write(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + '\t' + string + '\n')

class ProgressReporter(object):
    """
    Writes to out during iteration of unknown length
    """

    def __init__(self, name_of_process, report_step):

        self.name = name_of_process
        self.report_step = report_step
        self.position = int()
        self.next_limit = report_step

    def update(self):

        self.position += 1
        if self.position >= self.next_limit:
            report_progress(self.name + '\t' + "{:,}".format(self.position))
            self.next_limit += self.report_step

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

class FileReader(object):
    """
    Reads input files, handles gzip.
    """
    def __init__(self, filehandle, filehandle2=None):

        # Init variables setting
        self.filehandle = filehandle
        self.gzip = bool()

        # Open files as zipped or not not (depending on if they end with .gz)
        if self.filehandle[-3:] == '.gz':
            report_progress('File detected as gzipped, unzipping when reading')
            import gzip
            self.openfile = gzip.open(self.filehandle, 'r')
            self.gzip = True
        else:
            self.openfile = open(self.filehandle, 'r')

        # Paired end preparation
        self.filehandle2 = filehandle2
        if self.filehandle2:

            # Open files as zipped or not not (depending on if they end with .gz)
            if self.filehandle2[-3:] == '.gz':
                report_progress('File detected as gzipped, unzipping when reading')
                import gzip
                self.openfile2 = gzip.open(self.filehandle2, 'r')
            else:
                self.openfile2 = open(self.filehandle2, 'r')

    def fileReader(self):
        """
        Reads non-specific files as generator
        :return: lines
        """
        for line in self.openfile:
            if self.gzip:
                line = line.decode("utf-8")
            yield line

    def fastqReader(self):
        """
        Reads lines 4 at the time as generator
        :return: read as fastq object
        """

        line_chunk = list()
        for line in self.openfile:
            if self.gzip:
                line = line.decode("utf-8")
            line_chunk.append(line)
            if len(line_chunk) == 4:
                read = FastqRead(line_chunk)
                line_chunk = list()
                yield read

    def fastqPairedReader(self):
        """
        Reads two paired fastq files and returns a pair of two reads
        :return: read1 read2 as fastq read objects
        """

        line_chunk1 = list()
        line_chunk2 = list()
        for line1, line2 in zip(self.openfile, self.openfile2):
            if self.gzip:
                line1 = line1.decode("utf-8")
                line2 = line2.decode("utf-8")
            line_chunk1.append(line1)
            line_chunk2.append(line2)
            if len(line_chunk1) == 4 and len(line_chunk2) == 4:
                read1 = FastqRead(line_chunk1)
                read2 = FastqRead(line_chunk2)

                # Error handling
                if not read1.header.split()[0] == read2.header.split()[0]:
                    import sys
                    sys.exit('INPUT ERROR: Paired reads headers does not match.\nINPUT ERROR: Read pair number:\t'+str(progress.position+1)+'\nINPUT ERROR: '+str(read1.header)+'\nINPUT ERROR: '+str(read2.header)+'\nINPUT ERROR: Exiting')
                line_chunk1 = list()
                line_chunk2 = list()
                yield read1, read2

    def close(self):
        """
        Closes files properly so they can be re-read if need be.
        :return:
        """
        self.openfile.close()
        if self.filehandle2:
            self.openfile2.close()

class FastqRead(object):
    """
    Stores read as object.
    """

    def __init__(self, fastq_as_line):

        self.header = fastq_as_line[0].strip()
        self.seq = fastq_as_line[1].strip()
        self.comment = fastq_as_line[2].strip()
        self.qual = fastq_as_line[3].strip()

    def fastq_string(self):
        return self.header + '\n' + self.seq  + '\n' + self.comment  + '\n' + self.qual + '\n'

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
        import argparse
        global args

        parser = argparse.ArgumentParser(description="Tags bam files with barcode clustering information. Looks for raw "
                                                     "sequence in read header and puts barcode cluster ID in BC tag as "
                                                     "well as in header.")

        # Arguments
        parser.add_argument("input_mapped_bam", help=".bam file with mapped reads which is to be tagged with barcode id:s.")
        parser.add_argument("input_clstr", help=".clstr file from cdhit clustering.")
        parser.add_argument("output_tagged_bam", help=".bam file with barcode cluster id in the bc tag.")

        # Options
        parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                           "Not recommended due to different function "
                                                                           "names in python 2 and 3.")
        parser.add_argument("-s", "--skip_nonclust", action="store_true", help="Does not give cluster ID:s to clusters "
                                                                               "made out by only one sequence.")

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
