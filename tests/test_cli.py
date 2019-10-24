from pathlib import Path

from blr.cli.init import init


def test_init(tmpdir):
    init(tmpdir / "analysis", Path("testdata/reads.1.fastq.gz"))
