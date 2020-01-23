"""
Filter BAM for select SAM tag values.
"""
import logging
from collections import Counter
from tqdm import tqdm
import os

from blr.utils import get_bamtag, print_stats, PySAMIO

logger = logging.getLogger(__name__)


def main(args):
    logger.info("Starting")

    summary = Counter()
    values = get_values(args.values)
    counts = {value: 0 for value in values}
    summary["Values to collect"] = len(values)

    with PySAMIO(args.input, args.output, __name__) as (openin, openout):
        for read in tqdm(openin.fetch(until_eof=True), desc="Processing reads"):
            summary["Reads in"] += 1
            value = get_bamtag(read, args.tag)

            if value in values:
                summary["Reads out"] += 1
                counts[value] += 1
                openout.write(read)

    print_counts(counts)
    print_stats(summary, __name__)

    logger.info("Finished")


def get_values(values):
    if len(values) == 1 and os.path.isfile(values[0]):
        logger.info(f"Collecting values from file: {values[0]}.")
        with open(values[0], 'r') as values_file:
            return {v.strip() for v in values_file if v != ""}
    else:
        return set(values)


def print_counts(counts):
    sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    width_value = max([len(v[0]) for v in sorted_counts])
    print(f"{'Value':<{width_value + 1}}{'Count':>12}")
    for value, count in sorted_counts:
        print(f"{value:<{width_value+1}}{count:>12,}")


def add_arguments(parser):
    parser.add_argument("input",
                        help="SAM/BAM file")
    parser.add_argument("tag",
                        help="SAM tag")
    parser.add_argument("values", nargs="+",
                        help="Values to match for output. Could be file with entries on separate lines.")

    parser.add_argument("-o", "--output", default="-",
                        help="Write output BAM to file rather then stdout.")
