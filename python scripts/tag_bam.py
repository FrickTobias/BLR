#! /usr/bin/env python2

def main():
    """Takes a fastq file barcode sequences in the header and writes a barcode fasta file with only unique entries. """

    #
    # Imports & globals
    #
    import multiprocessing, pysam
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
    # Data processing & writing output
    #

    # Generate read:cluster dictionary from concatenated .clstr file (stores in Summary instance)
    with open(args.input_clstr, 'r') as openInfile:
        readAndProcessClusters(openInfile)

    # Tag bam file and write to output file
    infile = pysam.AlignmentFile(args.input_mapped_bam, 'rb')
    out = pysam.AlignmentFile(args.output_tagged_bam, 'wb', template=infile)
    for read in infile.fetch(until_eof=True):
        read_bc = read.query_name.split('_')[-1]
        try: bc_id = summaryInstance.read_to_barcode_dict[read_bc]
        except KeyError:
            Summary.writeLog('KeyError: ' + str(read_bc))
        read.set_tag('RG', str(bc_id), value_type='Z') # Stores as string, makes duplicate removal possible. Can do it as integer as well.
        out.write(read)
    infile.close()
    out.close()

def readAndProcessClusters(openInfile):
    """ Reads clstr file and builds read:clusterId dict in Summary instance."""

    # Set clusterInstance for first loop
    for first_line in openInfile:
        clusterInstance = ClusterObject(clusterId=first_line)
        break

    for line in openInfile:

        # Reports cluster to master dict and start new cluster instance
        if line.startswith('>'):
            summaryInstance.updateReadToClusterDict(clusterInstance.barcode_to_bc_dict)
            clusterInstance = ClusterObject(clusterId=line)

        # Add accession entry for current cluster id
        else:
            clusterInstance.addRead(line)

    # Add last cluster to master dict
    summaryInstance.updateReadToClusterDict(clusterInstance.barcode_to_bc_dict)

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
        """ Main funcion for overview of what is run. """

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
        parser.add_argument("input_mapped_bam", help=".bam file with mapped reads which is to be tagged with barcode id:s.")
        parser.add_argument("input_clstr", help=".clstr file from cdhit clustering.")
        parser.add_argument("output_tagged_bam", help=".bam file with barcode cluster id in the bc tag.")

        # Options
        parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                           "Not recommended due to different function "
                                                                           "names in python 2 and 3.")
        parser.add_argument("-p", "--processors", type=int, default=multiprocessing.cpu_count(),
                            help="Thread analysis in p number of processors. Example: python "
                                 "TagGD_prep.py -p 2 insert_r1.fq unique.fa")

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
    """ Summarizes chunks"""

    def __init__(self):
        self.read_to_barcode_dict = dict()
        self.CurrentClusterId = 0
        self.barcodeLength = int()
        log = args.output_tagged_bam.split('.')[:-1]
        self.log = '.'.join(log) + '.log'
        with open(self.log, 'w') as openout:
            pass

    def updateReadToClusterDict(self, input_dict):
        """ Merges cluster specific dictionaries to a master dictionary."""

        self.CurrentClusterId += 1
        for barcode in input_dict.keys():
            self.read_to_barcode_dict[barcode] = self.CurrentClusterId

    def writeLog(self, line):

        import time
        with open(self.log, 'a') as openout:
            openout.write(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + '\n')
            openout.write(line + '\n')

if __name__=="__main__": main()