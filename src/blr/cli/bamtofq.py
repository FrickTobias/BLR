"Writes fq files from bam files"

import logging
import pysam

from tqdm import tqdm
import dnaio

from blr import utils

logger = logging.getLogger(__name__)


def main(args):
    logger.info(f"Starting analysis")
    query_cache = dict()
    summary = Summary()
    out_interleaved = not args.r2_out
    with pysam.AlignmentFile(args.bam, "rb") as reader, \
            dnaio.open(args.r1_out, file2=args.r2_out, interleaved=out_interleaved, mode="w",
                       fileformat="fastq") as writer:
        for query in tqdm(reader):
            summary.tot_reads_in += 1

            # Read pair cache system
            if query.query_name not in query_cache:
                query_cache[query.query_name] = query
                continue
            else:
                query_1 = query_cache[query.query_name]
                query_2 = query
                del query_cache[query.query_name]

            # Fetch and add bam tag to header (if not found header is not modified)
            if args.tag:
                add_tag_to_header(query_1, query_2, args.tag)

            # Find out which query is read 1 and read 2 (so they will be output to correct file)
            r1, r2 = turn_readpair(query_1, query_2)

            r1_as_dnaio_object = dnaio._core.Sequence(name=r1.query_name, sequence=r1.seq, qualities=r1.qual)
            r2_as_dnaio_object = dnaio._core.Sequence(name=r2.query_name, sequence=r2.seq, qualities=r2.qual)

            writer.write(r1_as_dnaio_object, r2_as_dnaio_object)
            summary.tot_reads_out += 2

    summary.print_stats()
    logger.info(f"Finished")


def add_tag_to_header(query_1, query_2, tag):
    """
    Fetches a BAM tag and if present, appends to header of read.
    :param query_1: paired pysam alignment object
    :param query_2: paired pysam alignment object
    :return: None, objects are modified and not explicitly returned.
    """
    tag_string = utils.get_bamtag(query_1, tag)
    if tag_string:
        query_1.query_name = f"{query_1.query_name.rstrip()} {tag}:{tag_string}"
        query_2.query_name = query_1.query_name


def turn_readpair(query_1, query_2):
    """
    Finds and set which BAM alignment is read 1 and read 2. Returns None if not read 1 and read 2 found.
    :param query_1: paired pysam alignment object
    :param query_2: paired pysam alignment object
    :return: read1 and read2 as pysam alignment bjects
    """
    if query_1.is_read1 and query_2.is_read2:
        read1, read2 = query_1, query_2
    elif query_1.is_read2 and query_2.is_read1:
        read1, read2 = query_2, query_1
    else:
        read1, read2 = None, None

    return read1, read2


class Summary:

    def __init__(self):
        self.tot_reads_in = int()
        self.tot_reads_out = int()

    def print_stats(self):
        for object_variable, value in vars(self).items():
            logger.info(f"{object_variable}: {value}")


def add_arguments(parser):
    parser.add_argument("bam", help="BAM file input")
    parser.add_argument("r1_out", help="Name for r1 .fq output")
    parser.add_argument("r2_out", help="Name for r2 .fq output")
    parser.add_argument("-t", "--tag", type=str, metavar="<STRING>", help="BAM tag to add to header in fq out.")