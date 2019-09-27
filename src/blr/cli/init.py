"""
Create and initialize a new analysis directory.
"""
import logging
import os
import os.path
import sys
from pathlib import Path
from importlib_resources import read_binary

from ..utils import guess_paired_path

logger = logging.getLogger(__name__)


CONFIGURATION_FILE_NAME = "config.yaml"


def add_arguments(parser):
    parser.add_argument(
        "--reads1",
        "--r1",
        type=Path,
        required=True,
        metavar="READS",
        help="First paired-end read file (.fastq.gz). The second is found automatically.",
    )
    parser.add_argument("directory", type=Path, help="New analysis directory to create")


def main(args):
    init(args.directory, args.reads1)


def init(directory: Path, reads1: Path):
    if " " in str(directory):
        logger.error("The name of the analysis directory must not contain spaces")
        sys.exit(1)

    if directory.exists():
        logger.error(f"The target directory {directory!r} already exists.")
        sys.exit(1)

    fail_if_inaccessible(reads1)
    reads2 = guess_paired_path(reads1)
    if reads2 is None:
        logger.error("Could not determine second file of paired-end reads")
        sys.exit(1)
    fail_if_inaccessible(reads2)

    create_and_populate_analysis_directory(directory, reads1, reads2)

    logger.info(f"Directory {directory} initialized.")
    logger.info(
        'Edit %s/%s, then run "cd %s && blr run" to start the analysis',
        directory,
        CONFIGURATION_FILE_NAME,
        directory,
    )


def create_and_populate_analysis_directory(directory: Path, reads1: Path, reads2: Path):
    try:
        directory.mkdir()
    except OSError as e:
        logger.error(e)
        sys.exit(1)

    # Write the configuration file
    configuration = read_binary("blr", CONFIGURATION_FILE_NAME)
    with (directory / CONFIGURATION_FILE_NAME).open("wb") as f:
        f.write(configuration)

    create_symlink(reads1, directory, "reads.1.fastq.gz")
    create_symlink(reads2, directory, "reads.2.fastq.gz")


def fail_if_inaccessible(path):
    try:
        with path.open():
            pass
    except OSError as e:
        logger.error("Could not open %r: %s", path, e)
        sys.exit(1)


def create_symlink(readspath, dirname, target):
    if not os.path.isabs(readspath):
        src = os.path.relpath(readspath, dirname)
    else:
        src = readspath
    os.symlink(src, os.path.join(dirname, target))
