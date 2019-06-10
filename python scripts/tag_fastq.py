#! /usr/bin/env python

def main():
    """ Takes a .fastq file with trimmed reads and a .fasta file with clustered barcodes. Output is a tagged .fastq file. """

    global args, summaryObject
    from Bio import SeqIO
    import logging

    configureLogging('info')

    argParser() # Argument parsing
    logging.info('Arguments read successfully')

    summaryObject = Sum()

    readClusters(args.input_clstr)

    logging.info('Cluster file processed successfully')

    with open(args.input_trimmed_fastq, 'r') as openin, open(args.output_tagged_fastq, 'w') as openout:

        for record in SeqIO.parse(openin,'fastq'):

            name,bc_seq = record.id.split('_')

            try:
                consensus_seq = summaryObject.master_barcode_dict[bc_seq]
            except KeyError:
                continue

            record.id = name + '_' + consensus_seq
            record.description = record.id + ' BC:Z:' + consensus_seq + '-1'

            SeqIO.write(record,openout,'fastq')

    logging.info('Tagging completed')

class argParser(object):
    """Reads input arguments"""

    def __init__(self):

        argParser.parse(self)

    def parse(self):

        import argparse
        global args

        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("input_trimmed_fastq", help=".fastq file with trimmed reads which is to be tagged with barcode id:s.")
        parser.add_argument("input_clstr", help=".clstr file from cdhit clustering.")
        parser.add_argument("output_tagged_fastq", help=".fastq file tagged with barcode consensus sequence")
        args = parser.parse_args()

class configureLogging(object):

    import logging
    global logging

    def __init__(self,logType):

        log_format = ' %(asctime)s %(levelname)s: < %(message)s >'

        if logType == 'info':
            configureLogging.INFO(self,log_format)
        elif logType == 'warning':
            configureLogging.WARNING(self,log_format)
        elif logType == 'error':
            configureLogging.ERROR(self,log_format)
        elif logType == 'critical':
            configureLogging.CRITICAL(self,log_format)
        elif logType == 'debug':
            configureLogging.DEBUG(self,log_format)

    def INFO(self,log_format):

        logging.basicConfig(level=logging.INFO, format=log_format)

    def WARNING(self,log_format):

        logging.basicConfig(level=logging.WARNING, format=log_format)

    def ERROR(self,log_format):

        logging.basicConfig(level=logging.ERROR, format=log_format)

    def CRITICAL(self,log_format):

        logging.basicConfig(level=logging.CRITICAL, format=log_format)

    def DEBUG(self,log_format):

        logging.basicConfig(level=logging.DEBUG, format=log_format)


def readClusters(clstrfile):
    """ Reads clstr file and builds consensus_bc:bc dict in Sum"""

    new_cluster = False

    with open(clstrfile, 'r') as f:

        cluster = ClusterObject(clusterId=f.readline())  # set clusterId for the first cluster

        for line in f:
            # Reports cluster to master dict and start new cluster instance

            if line.startswith('>'):

                new_cluster = True

            elif new_cluster == True:

                new_cluster = False

                cluster = ClusterObject(clusterId=line)
                cluster.addBarcodeToDict(line)
                summaryObject.updateDict(cluster)

            else:
                cluster.addBarcodeToDict(line)

        # Add last cluster to master dict
        summaryObject.updateDict(cluster)


class ClusterObject(object):
    """ Cluster object, one for each cluster"""

    def __init__(self, clusterId):
        self.consensus_to_bc_dict = dict()
        self.consensus = clusterId.split()[-2].split(':')[-1].rstrip('.')

    def addBarcodeToDict(self, line):
        '''Add non-consensus bc to consensus_to_bc_dict'''

        accession = line.split()[2].rstrip('.')
        barcode = accession.split(':')[-1]
        self.consensus_to_bc_dict[barcode] = self.consensus

class Sum(object):

    def __init__(self):
        self.master_barcode_dict = dict()

    def updateDict(self, input_object):
        """ Merges cluster-specific dictionaries to a master dictionary."""

        input_dict = input_object.consensus_to_bc_dict

        for barcode in input_dict.keys():
            self.master_barcode_dict[barcode] = input_object.consensus

if __name__=="__main__": main()