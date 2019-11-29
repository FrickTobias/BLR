from pathlib import Path
import pysam
import pytest
import dnaio

from blr.cli.init import init
from blr.cli.run import run
from blr.cli.config import change_config

TESTDATA_READS = Path("testdata/reads.1.fastq.gz")
DEFAULT_CONFIG = "blr.yaml"
REFERENCE_GENOME = str(Path("testdata/chr1mini.fasta").absolute())
REFERENCE_VARIANTS = str(Path("testdata/HG002_GRCh38_GIAB_highconf.chr1mini.vcf").absolute())


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


def test_init(tmpdir):
    init(tmpdir / "analysis", TESTDATA_READS)


def test_config(tmpdir):
    workdir = tmpdir / "analysis"
    init(workdir, TESTDATA_READS)
    change_config(workdir / "blr.yaml", [("read_mapper", "bwa")])


@pytest.mark.parametrize("read_mapper", ["bwa", "bowtie2", "minimap2"])
def test_mappers(tmpdir, read_mapper):
    workdir = tmpdir / "analysis"
    init(workdir, TESTDATA_READS)
    change_config(
        workdir / DEFAULT_CONFIG,
        [("genome_reference", REFERENCE_GENOME), ("read_mapper", read_mapper)]
    )
    run(workdir=workdir, targets=["mapped.sorted.bam"])
    n_input_fastq_reads = 2 * count_fastq_reads(Path(workdir / "trimmed_barcoded.1.fastq.gz"))
    assert n_input_fastq_reads <= count_bam_alignments(workdir / "mapped.sorted.bam")


@pytest.mark.parametrize("duplicate_marker", ["sambamba", "samblaster"])
def test_duplicate_markers(tmpdir, duplicate_marker):
    workdir = tmpdir / "analysis"
    init(workdir, TESTDATA_READS)
    change_config(
        workdir / DEFAULT_CONFIG,
        [("genome_reference", REFERENCE_GENOME), ("duplicate_marker", duplicate_marker)]
    )
    run(workdir=workdir, targets=["mapped.sorted.tag.mkdup.bam"])
    n_input_fastq_reads = 2 * count_fastq_reads(Path(workdir / "trimmed_barcoded.1.fastq.gz"))
    assert n_input_fastq_reads <= count_bam_alignments(workdir / "mapped.sorted.tag.mkdup.bam")


def test_final_compressed_reads_exist(tmpdir):
    workdir = tmpdir / "analysis"
    init(workdir, TESTDATA_READS)
    change_config(
        workdir / DEFAULT_CONFIG,
        [("genome_reference", REFERENCE_GENOME)]
    )
    targets = ("reads.1.final.fastq.gz", "reads.2.final.fastq.gz")
    run(workdir=workdir, targets=targets)
    for filename in targets:
        assert Path(workdir / filename).exists()


def test_link_reference_variants(tmpdir):
    workdir = tmpdir / "analysis"
    init(workdir, TESTDATA_READS)
    change_config(
        workdir / DEFAULT_CONFIG,
        [("genome_reference", REFERENCE_GENOME), ("reference_variants", REFERENCE_VARIANTS)]
    )
    target = "variants.reference.vcf"
    run(workdir=workdir, targets=[target])
    assert Path(workdir / target).is_symlink()


def test_call_variants(tmpdir):
    workdir = tmpdir / "analysis"
    init(workdir, TESTDATA_READS)
    change_config(
        workdir / DEFAULT_CONFIG,
        [("genome_reference", REFERENCE_GENOME), ("reference_variants", "null")]
    )
    target = "variants.called.vcf"
    run(workdir=workdir, targets=[target])
    assert Path(workdir / target).is_file()
