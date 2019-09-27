import re
from pathlib import Path


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
