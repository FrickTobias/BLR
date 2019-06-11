#! /usr/bin/env python2

def main():
    """Takes a fastq file barcode sequences in the header and writes a barcode fasta file with only unique entries. """

    #
    # Imports & globals
    #
    global args, output_tagged_bamfile, sys, time, progress
    import sys, time

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
    cluster_dict = readAndProcessClusters(clstr_generator.openfile)
    clstr_generator.close()

    # Input clustering information in reads
    progress = ProgressReporter("Read pairs processed", 1000000)
    reads_written = int()
    reads_with_cluster_info = int()
    readpair_generator = FileReader(args.input_fastq1, args.input_fastq2)
    with open(args.output1, 'w') as openout1, open(args.output2, 'w') as openout2:
        for read1, read2 in readpair_generator.fastqPairedReader():
            # Fetch bc and if clustered, add barcode cluster ID to header
            bc_seq = read1.header.split()[0].split('_')[1]
            if bc_seq in cluster_dict:
                cluster_ID = cluster_dict[bc_seq]

                if not args.alternative_format:
                    read1.header = read1.header.split()[0] + ' BC:Z:' + str(cluster_ID) + '-1'
                    read2.header = read2.header.split()[0] + ' BC:Z:' + str(cluster_ID) + '-1'
                else:
                    read1.header = read1.header.split()[0] + '_BC:Z:' + str(cluster_ID) + ' ' + read1.header.split()[1]
                    read2.header = read2.header.split()[0] + '_BC:Z:' + str(cluster_ID) + ' ' + read2.header.split()[1]
                reads_with_cluster_info += 1

            # Write & report
            reads_written += 1
            openout1.write(read1.fastq_string())
            openout2.write(read2.fastq_string())
            progress.update()

    # Close input fastqs
    readpair_generator.close()

    # Reporting
    report_progress("Reads written to output file:\t" + str(reads_written))
    report_progress("Reads with clustering information:\t" + str(reads_with_cluster_info))
    report_progress("Finished")

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
                del cluster_dict[bc_seq]
            seqs_in_cluster = 0

            # Save consensus bc seq instead of barcode cluster ID
            if args.consensus_seq:
                cluster_ID = False
            # Runs an iterator giving a cluster ID to mark it is a cluster
            else:
                cluster_ID += 1
        else:
            # Fetch bc seq from cdhit output format
            bc_seq = line.split()[2].lstrip('>').rstrip('...').split(':')[2]
            seqs_in_cluster += 1

            # If not using cluster ID, take first seq (cluster consensus sequence) and use as value
            if not cluster_ID:
                cluster_ID = bc_seq
            cluster_dict[bc_seq] = cluster_ID

    if args.skip_nonclust and seqs_in_cluster < 2:
        del cluster_dict[bc_seq]

    return(cluster_dict)

def report_progress(string):
    """
    Writes a time stamp followed by a message (=string) to standard out.
    Input: String
    Output: [date]  string
    """
    sys.stderr.write(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + '\t' + string + '\n')

class ClusterObject(object):
    """
    Cluster object
    """

    currentClusterId = int()

    def __init__(self, clusterId):
        self.barcode_to_bc_dict = dict()
        ClusterObject.currentClusterId += 1
        self.Id = ClusterObject.currentClusterId  # Remove 'Cluster' string and \n from end
        self.numberOfReads = 0

    def addRead(self, line):
        accession = line.split()[2].rstrip('.')
        barcode = accession.split(':')[-1]
        self.barcode_to_bc_dict[barcode] = self.Id  # Extract header and remove '...'
        self.numberOfReads += 1
        if self.numberOfReads == 1:
            summaryInstance.updateClusterIdToBarcodeSeqDict(barcode, self.Id)


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

        parser = argparse.ArgumentParser(description="Tags a fastq files (with _BDVH...VH in header) with CD-HIT-454 "
                                                     "output (NNN.clstr file). Stores ist as ")

        # Arguments
        parser.add_argument("input_fastq1", help="r1 .fastq file with trimmed reads which is to be tagged with barcode id:s.")
        parser.add_argument("input_fastq2", help="r2 .fastq file with trimmed reads which is to be tagged with barcode id:s.")
        parser.add_argument("input_clstr", help=".clstr file from cdhit clustering.")
        parser.add_argument("output1", help="r1 .fastq file with barcode cluster id in the bc tag.")
        parser.add_argument("output2", help="r2 .fastq file with barcode cluster id in the bc tag.")

        # Options
        parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                           "Not recommended due to different function "
                                                                           "names in python 2 and 3.")
        parser.add_argument("-s", "--skip_nonclust", action="store_true", help="Does not give cluster ID:s to clusters "
                                                                                "made out by only one sequence.")
        parser.add_argument("-af", "--alternative_format", action="store_true", help="Writes in alternative format. Normal"
                                                                                   "format writes yields "
                                                                                   "'@header_BcSeq BC:Z:integer' (which "
                                                                                   "is compatible with BWA) and af "
                                                                                   "format yields "
                                                                                   "'@header_BcSeq_BC:Z:integer more_info'")
        parser.add_argument('-cs', '--consensus_seq', action="store_true", help="Writes consensus sequence of barcode "
                                                                                "instead of the barcode cluster ID.")

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