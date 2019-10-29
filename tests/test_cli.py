from pathlib import Path
import pysam
import pytest
import dnaio

from blr.cli.init import init
from blr.cli.run import run

TESTDATA_READS = Path("testdata/reads.1.fastq.gz")


def count_bam_alignments(path):
    with pysam.AlignmentFile(path) as af:
        n = 0
        for _ in af:
            n += 1
    return n


def count_fastq_reads(path):
    with dnaio.open(path) as f:
        n = 0
        for _ in f:
            n += 1
    return n


def copy_config(source, target, genome_reference=None, read_mapper=None, duplicate_marker=None):
    """Copy config, possibly changing genome_reference or read_mapper"""

    with open(source) as infile:
        with open(target, "w") as outfile:
            for line in infile:
                if genome_reference is not None and line.startswith("genome_reference:"):
                    path = Path("testdata/chr1mini.fasta").absolute()
                    line = f"genome_reference: {path}\n"
                if read_mapper is not None and line.startswith("read_mapper:"):
                    line = f"read_mapper: {read_mapper}\n"
                if duplicate_marker is not None and line.startswith("duplicate_marker:"):
                    line = f"duplicate_marker: {duplicate_marker}\n"
                outfile.write(line)


def test_init(tmpdir):
    init(tmpdir / "analysis", TESTDATA_READS)


@pytest.mark.parametrize("read_mapper", ["bwa", "bowtie2", "minimap2"])
def test_mappers(tmpdir, read_mapper):
    workdir = tmpdir / "analysis"
    init(workdir, TESTDATA_READS)
    copy_config(
        "tests/test_config.yaml",
        workdir / "blr.yaml",
        genome_reference=str(Path("testdata/chr1mini.fasta").absolute()),
        read_mapper=read_mapper,
    )
    run(workdir=workdir, targets=["mapped.sorted.bam"])
    n_input_fastq_reads = 2 * count_fastq_reads(Path(workdir / "trimmed_barcoded.1.fastq.gz"))
    assert n_input_fastq_reads <= count_bam_alignments(workdir / "mapped.sorted.bam")


@pytest.mark.parametrize("duplicate_marker", ["sambamba", "samblaster"])
def test_duplicate_markers(tmpdir, duplicate_marker):
    workdir = tmpdir / "analysis"
    init(workdir, TESTDATA_READS)
    copy_config(
        "tests/test_config.yaml",
        workdir / "blr.yaml",
        genome_reference=str(Path("testdata/chr1mini.fasta").absolute()),
        read_mapper="bwa",
        duplicate_marker=duplicate_marker
    )
    run(workdir=workdir, targets=["mapped.sorted.tag.mkdup.bam"])
    n_input_fastq_reads = 2 * count_fastq_reads(Path(workdir / "trimmed_barcoded.1.fastq.gz"))
    assert n_input_fastq_reads <= count_bam_alignments(workdir / "mapped.sorted.tag.mkdup.bam")
