"""
Removes barcode duplicates (two different barcode sequences origin to the same droplet, tagging the same tagmented long
molecule) by merging barcode sequences for reads sharing duplicates.

Condition to call barcode duplicate:

Two positions (positions defined as a unique set of read_start, read_stop, mate_start, mate_stop))
at a maximum of W (--window, default 100kbp, between = max(downstream_pos)-max(downstream_pos)) bp apart
sharing more than one barcode (share = union(bc_set_pos1, bc_set_pos2)).
"""

import pysam
import logging
from tqdm import tqdm

logger = logging.getLogger(__name__)


def main(args):
    logger.info("Starting Analysis")

    summary = Summary()

    current_cache_rp = dict()
    cache_reads = dict()
    chrom_prev = None
    pos_prev = None
    merge_dict = dict()
    cache_dup_pos = dict()
    with pysam.AlignmentFile(args.input_tagged_bam, "rb") as openin:
        for read in tqdm(openin.fetch(until_eof=True), desc="Processing reads"):
            summary.tot_reads += 1

            # Wait for mate of read until process
            if read.query_name not in cache_reads:
                cache_reads[read.query_name] = read
                continue
            else:
                mate = cache_reads[read.query_name]
                del cache_reads[read.query_name]

            # Requirements: read mapped, mate mapped and read has barcode tag
            rp_meet_requirements, bc_new, summary = meet_requirements(read=read, mate=mate, summary=summary,
                                                                      barcode_tag=args.barcode_tag)
            if rp_meet_requirements:
                chrom_new = read.reference_name
                summary.reads_analyzed += 2

                pos_new = read.reference_start
                rp_pos_tuple = (mate.reference_start, mate.reference_end, read.reference_start, read.reference_end)

                # Save all rp until position is new to know if any are marked as duplicates.
                if chrom_new == chrom_prev and pos_new == pos_prev:

                    if rp_pos_tuple in current_cache_rp:
                        current_cache_rp[rp_pos_tuple].add_read_pair_and_bc(read=read, mate=mate, bc=bc_new)
                    else:
                        current_cache_rp[rp_pos_tuple] = CacheReadPairTracker(rp_pos_tuple=rp_pos_tuple,
                                                                              chromosome=chrom_new, read=read,
                                                                              mate=mate, bc=bc_new)

                # New pos tuple => Reset rp cache tracker
                else:
                    # If duplicates are found, try and find bc duplicates
                    for cache_read_pair_tracker in current_cache_rp.values():
                        if cache_read_pair_tracker.duplicate_read_pair():
                            summary.reads_at_analyzed_dup_position += len(cache_read_pair_tracker.current_reads) * 2
                            merge_dict, cache_dup_pos = seed_duplicates(
                                merge_dict=merge_dict,
                                cache_dup_pos=cache_dup_pos,
                                pos_new=cache_read_pair_tracker.position_tuple_ID,
                                bc_new=cache_read_pair_tracker.barcodes,
                                window=args.window)
                    current_cache_rp = dict()
                    current_cache_rp[rp_pos_tuple] = CacheReadPairTracker(rp_pos_tuple=rp_pos_tuple,
                                                                          chromosome=chrom_new, read=read, mate=mate,
                                                                          bc=bc_new)

                    # New chr => reset dup pos cache and rp cache tracker
                    if not chrom_new == chrom_prev:
                        cache_dup_pos = dict()

                pos_prev = pos_new
                chrom_prev = chrom_new

        summary.unmapped_reads += (len(cache_reads))

        # Last chunk
        for cache_read_pair_tracker in current_cache_rp.values():
            if cache_read_pair_tracker.duplicate_read_pair():
                summary.reads_at_analyzed_dup_position += len(cache_read_pair_tracker.current_reads) * 2
                merge_dict, cache_dup_pos = seed_duplicates(merge_dict=merge_dict, cache_dup_pos=cache_dup_pos,
                                                            pos_new=cache_read_pair_tracker.position_tuple_ID,
                                                            bc_new=cache_read_pair_tracker.barcodes,
                                                            window=args.window)

    # Remove several step redundancy (5 -> 3, 3 -> 1) => (5 -> 1, 3 -> 1)
    reduce_several_step_redundancy(merge_dict)
    summary.barcodes_removed = len(merge_dict)

    # Write outputs
    bc_seq_already_written = set()
    with open(args.merge_log, "w") as bc_merge_file, \
            pysam.AlignmentFile(args.input_tagged_bam, "rb") as infile, \
            pysam.AlignmentFile(args.output_bam, "wb", template=infile) as out:
        for read in tqdm(infile.fetch(until_eof=True), desc="Writing output", total=summary.tot_reads):

            # If read barcode in merge dict, change tag and header to compensate.
            try:
                previous_barcode_id = int(read.get_tag(args.barcode_tag))
            except KeyError:
                previous_barcode_id = None
            if previous_barcode_id in merge_dict:
                summary.reads_with_new_tag += 1
                new_barcode_id = str(merge_dict[previous_barcode_id])
                read.set_tag(args.barcode_tag, new_barcode_id, value_type="Z")
                read.query_name = "_".join(
                    read.query_name.split("_")[:-1]) + "_" + args.barcode_tag + ":Z:" + new_barcode_id

                # Merge file writing
                if new_barcode_id not in bc_seq_already_written:
                    bc_seq_already_written.add(new_barcode_id)
                    bc_merge_file.write(f"{previous_barcode_id},{new_barcode_id}\n")

            out.write(read)

    summary.print_stats(barcode_tag=args.barcode_tag)
    logger.info("Finished")


def meet_requirements(read, mate, summary, barcode_tag):
    """
    Checks so read pair meets requirements before being used in analysis.
    :param read: pysam read
    :param mate: pysam mate
    :return: bool
    """

    rp_meet_requirements = True
    bc_new = None

    if read.is_unmapped or mate.is_unmapped:
        if read.is_unmapped != mate.is_unmapped:
            summary.unmapped_mates += 1
            summary.unmapped_reads += 1
        else:
            summary.unmapped_reads += 2
        rp_meet_requirements = False

    try:
        bc_new = int(read.get_tag(barcode_tag))
    except KeyError:
        summary.non_tagged_reads += 2
        rp_meet_requirements = False

    return rp_meet_requirements, bc_new, summary


class CacheReadPairTracker:
    """
    Stores read pairs and keeps track is reads/mates are marked as duplicate for that set of reads.
    """

    def __init__(self, rp_pos_tuple, chromosome, read, mate, bc):

        self.position_tuple_ID = rp_pos_tuple
        self.chromosome = chromosome
        self.current_reads = dict()
        self.barcodes = set()
        self.read_pos_has_duplicates = bool()
        self.mate_pos_has_duplciates = bool()

        self.add_read_pair_and_bc(read=read, mate=mate, bc=bc)

    def add_read_pair_and_bc(self, read, mate, bc):

        if read.is_duplicate:
            self.read_pos_has_duplicates = True
        if mate.is_duplicate:
            self.mate_pos_has_duplciates = True

        self.current_reads[read.header] = (read, mate)
        self.barcodes.add(bc)

    def duplicate_read_pair(self):

        return self.read_pos_has_duplicates and self.mate_pos_has_duplciates


def seed_duplicates(merge_dict, cache_dup_pos, pos_new, bc_new, window):
    """
    Builds up a merge dictionary for which any keys should be overwritten by their value. Also keeps all previous
    positions saved in a cache dict ([pos_tuple] = bc_set) in which all reads which still are withing the window
    size are saved.
    :param merge_dict: Dict for tracking which bc_ids shuold be merged (directional) .
    :param cache_dup_pos: Dict for tracking previous duplicate positions and their barcode sets.
    :param pos_new: Position to be analyzed and subsequently saved to cache.
    :param bc_new: Set of BC:s at new pos
    :param window: Max distance allowed between prev read tuple stop and new read tuple start to call bc dup.
    :return: updated merge_dict & updated cache_dup_pos
    """

    if len(bc_new) >= 2:
        pos_start_new = min(pos_new)
        for pos_prev, bc_prev in cache_dup_pos.copy().items():
            pos_stop_prev = max(pos_prev)
            if pos_stop_prev + window >= pos_start_new:
                bc_union = bc_new & bc_prev

                # If two or more unique bc ID:s are found, add [big_clust_ID] = smallest_clust_ID to merge dict
                if len(bc_union) >= 2:

                    bc_union_sort = sorted(bc_union)
                    bc_union_min = bc_union_sort[0]
                    for bc_other in bc_union_sort[1:]:

                        # Never add give one key multiple values (connect the prev/new values instead)
                        if bc_other in merge_dict:
                            bc_min_prev = find_min_bc(bc_other, merge_dict)
                            bc_union_real_min = find_min_bc(bc_union_min, merge_dict)
                            if not bc_min_prev == bc_union_real_min:
                                if min(bc_union_real_min, bc_min_prev) == bc_union_real_min:
                                    merge_dict[bc_min_prev] = bc_union_real_min
                                else:
                                    merge_dict[bc_union_real_min] = bc_min_prev

                        # Normal case, just add high_bc_id => min_bc_id
                        else:
                            merge_dict[bc_other] = bc_union_min
            else:
                del cache_dup_pos[pos_prev]  # remove positions outside of window since pos are sorted
        cache_dup_pos[pos_new] = bc_new

    return merge_dict, cache_dup_pos


def find_min_bc(bc_minimum, merge_dict):
    """
    Goes through merge dict and finds the lowest current value for a chain of key-value entries. E.g if merge_dict has
    100 => 80, 80 => 60, 60 => 40 it will return 40 if any of the values 100,80,60 or 40 are given.
    :return: lowest clstr id for key-value chain
    """

    while True:
        try:
            bc_minimum = merge_dict[bc_minimum]
        except KeyError:
            break

    return bc_minimum


def reduce_several_step_redundancy(merge_dict):
    """
    Takes translation  dict saved in object and makes sure 5->3, 3->1 becomes 5->1, 3->1
    :param merge_dict: "messy" merge_dict
    :return: "clean" merge_dict
    """

    for barcode_to_remove in sorted(merge_dict.keys()):
        real_min = find_min_bc(barcode_to_remove, merge_dict)
        if not real_min == merge_dict[barcode_to_remove]:
            del merge_dict[barcode_to_remove]
            merge_dict[barcode_to_remove] = real_min

    return merge_dict


class Summary:

    def __init__(self):
        self.tot_reads = int()
        self.reads_with_new_tag = int()
        self.unmapped_reads = int()
        self.unpaired_reads = int()
        self.non_tagged_reads = int()
        self.reads_analyzed = int()
        self.unmapped_mates = int()
        self.reads_at_analyzed_dup_position = int()
        self.barcodes_removed = int()

    def print_stats(self, barcode_tag):
        """
        Prints stats to terminal
        """

        # Read stats
        logger.info(f"Total Reads in file: {self.tot_reads}")
        logger.info(f"  Unmapped reads: {self.unmapped_reads}")
        logger.info(f"  Reads with unmapped mate: {self.unmapped_mates}")
        logger.info(f"  Reads without {barcode_tag} tag: {self.non_tagged_reads}")
        logger.info(f"  Reads analyzed: {self.reads_analyzed}")
        logger.info(f"Reads at dup pos: {self.reads_at_analyzed_dup_position}")
        logger.info(f"Reads with new tag: {self.reads_with_new_tag}")
        logger.info(f"{barcode_tag} tags removed: {self.barcodes_removed}")


def add_arguments(parser):
    parser.add_argument("input_tagged_bam", help="Sorted .bam file tagged with barcodes (-bc).")
    parser.add_argument("output_bam", help="Sorted .bam file without barcode duplicates.")
    parser.add_argument(
        "merge_log",
        help=".csv log file containing all merges done. File is in format: {old barcode id},{new barcode id}")
    parser.add_argument(
        "-bc", "--barcode_tag", metavar="<BARCODE_TAG>", type=str, default="BC",
        help=".bam file tag in which the barcode is specified in. DEFAULT: BC")
    parser.add_argument(
        "-w", "--window", metavar="<INTEGER>", type=int, default=100000,
        help="Window size. Duplicate positions within this distance will be used to find cluster duplicates.")
