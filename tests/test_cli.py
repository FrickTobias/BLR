from pathlib import Path

from blr.cli.init import init
from blr.cli.run import run


def copy_config(source, target, genome_reference=None):
    """Copy config, possibly changing the genome_reference"""

    with open(source) as infile:
        with open(target, "w") as outfile:
            for line in infile:
                if genome_reference is not None and line.startswith("genome_reference:"):
                    line = "genome_reference: " + str(Path("testdata/chr1mini.fasta").absolute())
                outfile.write(line)


def test_init(tmpdir):
    init(tmpdir / "analysis", Path("testdata/reads.1.fastq.gz"))


def test_run(tmpdir):
    workdir = tmpdir / "analysis"
    init(workdir, Path("testdata/reads.1.fastq.gz"))
    copy_config(
        "tests/test_config.yaml",
        workdir / "blr.yaml",
        genome_reference=str(Path("testdata/chr1mini.fasta").absolute()),
    )
    run(workdir=workdir, targets=["mapped.sorted.bam"])
