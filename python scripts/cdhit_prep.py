#! /usr/bin python3

def main():
    """Takes a fastq file barcode sequences in the header and writes a barcode fasta file with only unique entries. """

    #
    # Imports & globals
    #
    import multiprocessing, argparse, sys
    global args

    #
    # Argument parsing
    #
    parser = argparse.ArgumentParser(description=__doc__)
    # Arguments
    parser.add_argument("input_fastq",help="Read file with barcode sequences as last element of accession row, separated "
                                           "by and underline. Example: '@ACCESSION_AGGTCGTCGATC'. Also handles "
                                           "'@ACCESSSION_AGGTCGTCGATC MORE_ACCESSION'.")
    parser.add_argument("output_fasta",help="Output file name with unique barcode sequences.")
    # Options
    parser.add_argument("-f","--filter", type=int, default=2, help="Filter file for minimum amount of reads default=2")
    parser.add_argument("-F","--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                      "Not recommended due to different function "
                                                                      "names in python 2 and 3.")
    parser.add_argument("-p","--processors", type=int, default=multiprocessing.cpu_count(), help="Thread analysis in p number of processors. Example: python "
                                                            "TagGD_prep.py -p 2 insert_r1.fq unique.fa")
    parser.add_argument("-r", "--reduce_complexity", type=int, help="Reduce complexity of clustering by using the "
                                                                    "first r bases as index for dividing barcodes into "
                                                                    "4^r number of files. Will then use output to arg "
                                                                    "as output directory instead of output file.")
    parser.add_argument("-s", "--space_separation", action="store_true", help="Alternative fastq input format, requires "
                             "barcode sequence to be separated by space "
                             "and written as last element of accession.")
    args = parser.parse_args()

    #
    # Version control
    #
    if sys.version_info.major == 3:
        pass
    else:
        sys.stderr.write('\nWARNING: you are running python ' + str(sys.version_info.major) + ', this script is written for python 3.')
        if not args.force_run:
            sys.stderr.write('\nAborting analysis.')
            sys.exit()
        else:
            sys.stderr.write('\nForcing run. This might yield inaccurate results.\n')

    #
    # Processors
    #
    processor_count = args.processors
    max_processor_count = multiprocessing.cpu_count()
    if processor_count == max_processor_count:
        pass
    elif processor_count > max_processor_count:
        sys.stderr.write('Computer does not have ' + str(processor_count) + ' processors, running with default (' + str(max_processor_count) + ')\n')
        processor_count = max_processor_count
    else:
        sys.stderr.write('Running with ' + str(processor_count) + ' processors.\n')

    #
    # Filtering
    #
    if not args.filter == 2:
        sys.stderr.write('Filtering barcodes with less than ' + str(args.filter) + ' reads\n')

    #
    # Initials
    #
    summaryInstance = Summary()

    #
    # Data processing
    #
    poolOfProcesses = multiprocessing.Pool(processor_count, maxtasksperchild=100000000)
    # Runs optimized version of software
    if not args.space_separation:
        parallelResults = poolOfProcesses.imap_unordered(forEveryRead, readInsertFiles(args.input_fastq, 4) ,chunksize=10000)
    # If options are used, runs function which checks which are active
    else:
        parallelResults = poolOfProcesses.imap_unordered(optionalForEveryRead, readInsertFiles(args.input_fastq, 4),chunksize=10000)
    for barcode_dict in parallelResults:
        Summary.mergeDicts(summaryInstance, barcode_dict)

    #
    # Complexity reduction
    #
    if args.reduce_complexity:
        reduceComplexity(summaryInstance)

    else:

        #
        # Write outfile
        #
        Summary.writeOutput(summaryInstance, args.output_fasta)

    #
    # Stdout
    #
    print ('Number of barcodes with less than ' + str(summaryInstance.minReadCount) + ' reads: ' + str(summaryInstance.filteredReadCount))
    print ('(Out of ' + str(summaryInstance.totalReadCount) + ')')

def readInsertFiles(infile, chunkSize):
    """ File reader, reads file chunk by chunk."""

    from itertools import islice
    with open(infile,'r') as openInfile:

        # Loop until EOF (chunk == None)
        while True:
            chunk = list(islice(openInfile,chunkSize))
            if not chunk:
                break
            yield chunk

def forEveryRead(chunk):
    """ Function run for every read pair found in insert files."""

    # Define variables
    barcode_dict = dict()

    # Dividing chunk into reads
    read_list = extractReadFromChunk(chunk)

    for read in read_list:

        # Create read object
        readInstance = ReadObject(read)

        # Check if barcode sequence has occurred before
        try: barcode_dict[readInstance.barcode] += 1
        except KeyError:
            barcode_dict[readInstance.barcode] = 1

    #import sys
    #sys.stderr.write('10k reads processed...')

    return barcode_dict

def extractReadFromChunk(chunk):
    """ Creates reads from chunks of input files."""

    read_list = list()
    for read in [chunk[x:x+4] for x in range(0, len(chunk), 4)]:
        read_list.append(read)

    return read_list

def optionalForEveryRead(chunk):
    """The same as forEveryRead but it will check for options when running which might lower performance."""

    # Define variables
    barcode_dict = dict()

    # Dividing chunk into reads
    read_list = extractReadFromChunk(chunk)

    for read in read_list:

        # Create read object
        readInstance = ReadObject(read)
        # Overwrites barcode sequence to fit space separation format if -s option is used.
        if args.space_separation:
            ReadObject.optionSpaceSeparatino(readInstance, read)

        # Check if barcode sequence has occurred before
        try: barcode_dict[readInstance.barcode] += 1
        except KeyError:
            barcode_dict[readInstance.barcode] = 1

    return barcode_dict

def reduceComplexity(summaryInstance):
    """ Uses r first bases as indexes and divides files accordingly."""

    #
    # Imports
    #
    import os, sys

    #
    # Generate dict with possibilities indexes
    #
    bases_list = ['A','T','C','G']
    index_dict = {'':dict()}
    reduction = args.reduce_complexity
    # Repeat for lenth of r
    for i in range(reduction):
        # Extend every key...
        for key in index_dict.copy().keys():
            # ... with one of every base
            for base in bases_list:
                index_dict[key+base] = dict()
            # Remove non-extended key
            del index_dict[key]

    #
    # Classify reads into indexes
    #
    for barcode, count in summaryInstance.barcode_dict.items():
        try: index_dict[barcode[:reduction]][barcode] = count
        except KeyError:
            summaryInstance.notATGCindex.append(barcode)
    #
    # Create output directory
    #
    output_directory = args.output_fasta
    try: os.mkdir(output_directory)
    except OSError:
        pass

    #
    # Write output fasta files named after their indexes.
    #
    for index in index_dict.keys():
        summaryInstance.barcode_dict = index_dict[index]
        summaryInstance.writeOutput(output_directory + '/' + index + '.fa')

    not_ATGC_index = len(summaryInstance.notATGCindex)
    if not_ATGC_index > 0:
        sys.stderr.write(str(not_ATGC_index) + ' indexing barcode sequences found with something else than ATCG.\n')
        sys.stderr.write('Writing not_ATCG_index.txt file with those indexes...\n')
        with open(args.output_fasta + '/not_ATCG_index.txt', 'w') as openout:
            for index in summaryInstance.notATGCindex:
                openout.write(index + '\n')

class ReadObject(object):
    """ Read objects created for passing into the filtering."""

    def __init__(self, read):
        self.barcode = read[0].split()[0].split('_')[-1]

    def optionSpaceSeparatino(self, read):
        """ For -s option the barcode is at another position so the barcode is redefined."""
        self.barcode = read[0].split()[-1].rstrip()

class Summary(object):
    """ Summarizes chunks"""

    def __init__(self):
        self.barcode_dict = dict()
        self.filteredReadCount = 0
        self.totalReadCount = 0
        self.minReadCount = args.filter
        self.notATGCindex = list()

    def mergeDicts(self, barcode_dict):

        for new_key in barcode_dict.keys():
            try: self.barcode_dict[new_key] += barcode_dict[new_key]
            except KeyError:
                self.barcode_dict[new_key] = barcode_dict[new_key]

    def writeOutput(self, output_fasta):
        """ Writes output file."""
        bc_id = 1
        with open(output_fasta, 'w') as output:
            for barcode in self.barcode_dict.keys():
                if self.barcode_dict[barcode] >= self.minReadCount:
                    output.write('>'+str(bc_id) + ':' + str(self.barcode_dict[barcode]) + ':' + barcode + '\n' + barcode +'\n')
                    bc_id += 1
                    self.totalReadCount += self.barcode_dict[barcode]
                else:
                    self.filteredReadCount += self.barcode_dict[barcode]
                    self.totalReadCount += self.barcode_dict[barcode]

if __name__ == "__main__": main()
