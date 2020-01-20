import re
from pathlib import Path
import sys
import pysam


def is_1_2(s, t):
    """
    Determine whether s and t are identical except for a single character of
    which one of them is '1' and the other is '2'.
    """
    differences = 0
    one_two = {"1", "2"}
    for c1, c2 in zip(s, t):
        if c1 != c2:
            differences += 1
            if differences == 2:
                return False
            if {c1, c2} != one_two:
                return False
    return differences == 1


def guess_paired_path(path: Path):
    """
    Given the path to a file that contains the sequences for the first read in a
    pair, return the file that contains the sequences for the second read in a
    pair. Both files must have identical names, except that the first must have
    a '1' in its name, and the second must have a '2' at the same position.

    Return None if no second file was found or if there are too many candidates.

    >>> guess_paired_path(Path('file.1.fastq.gz'))  # doctest: +SKIP
    'file.2.fastq.gz'  # if that file exists
    """
    name = path.name
    # All lone 1 digits replaced with '?'
    name_with_globs = re.sub(r"(?<![0-9])1(?![0-9])", "?", name)
    paths = [p for p in path.parent.glob(name_with_globs) if is_1_2(str(p), str(path))]
    if len(paths) == 1:
        return paths[0]
    return None


def get_bamtag(pysam_read, tag):
    """
    Fetches tags from bam files. Return an empty value of the same type if not found.
    :param pysam_read: pysam read object
    :param tag: bam tag to fetch
    :return: bam tag value
    """
    try:
        return pysam_read.get_tag(tag)
    except KeyError:
        return None


def print_stats(summary, name=None, value_width=15, print_to=sys.stderr):
    """
    Prints stats in nice table with two column for the key and value pairs in summary
    :param summary: collections.Coutner object
    :param name: name of script for header e.g. '__name__'
    :param value_width: width for values column in table
    :param print_to: Where to direct output
    """
    # Get widths for formatting
    max_name_width = max(map(len, summary.keys()))
    width = value_width + max_name_width + 1

    # Header
    print("="*width, file=print_to)
    print(f"STATS SUMMARY - {name}", file=print_to)
    print("-"*width, file=print_to)

    # Print stats in columns
    for name, value in summary.items():
        print(f"{name:<{max_name_width}} {value:>{value_width},}", file=print_to)

    print("="*width, file=print_to)


def create_header(file, name):
    """
    Create SAM header dict with new tool and command line argument information based on template file. Appends new PG
    entry with tool name (ID), software name (PN), command line arguments (CL) to track the tools applied to the file.
    Use in output SAM/BAM file as 'header' attribute.

    :param file: string. Path to input SAM/BAM file.
    :param name: string. Pass '__name__' variable to be used to get program and tool name.
    :return: pysam.AlignmentHeader object
    """
    def check_name(identifier, prev_entries):
        """
        Check that identifier does not appear in previous header entries. Update with numbered suffix if the name is
        not unique.
        :param identifier: string. Program record identifier
        :param prev_entries: dict. Dictionary of "PG" header entries from pysam.AlignmentFile header.
        :return: string. Updated identifier
        """
        nr = 0
        updated_identifier = identifier
        while any(updated_identifier == e["ID"] for e in prev_entries):
            nr += 1
            updated_identifier = "_".join([identifier, str(nr)])
        return updated_identifier

    id_name = name.split(".")[-1]
    program_name = name.split(".")[0]
    cmd_line = f"\"{' '.join(sys.argv)}\""

    with pysam.AlignmentFile(file, "rb") as template:
        header = template.header.to_dict()
        pg_entries = header["PG"]

        # Check that id is unique and does not appear in header, otherwise add number at end.
        id_name = check_name(id_name, pg_entries)

        # TODO add version information (VN tag).
        pg_entries.append({
            "ID": id_name,          # Program record identifier. Must be unique
            "PN": program_name,     # Program name
            "CL": cmd_line          # Command line
        })

        header["PG"] = pg_entries

    return pysam.AlignmentHeader.from_dict(header)
