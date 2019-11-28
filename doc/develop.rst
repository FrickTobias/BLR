Development
===========

Conda environment files
-----------------------

There are two types of files that describe Conda environments.

- The file ``environment.yml`` contains abstract dependencies such as ``pysam`` or
  ``bowtie2``. This file is managed manually and needs to be
  updated whenever there are new dependencies or when the required version for a
  dependency changes.

- The ``environment.linux.lock.yml`` and ``environment.osx.lock.yml`` files
  (lock files) contain a fully specified description of the entire environment,
  with locked-down versions.  These files are used to create the test
  environment.

Use the script ``misc/condalock.sh`` to update the lock files whenever you make
changes to ``environment.yml``.


SAM Tags
--------
Specifications on SAM-tags used for holding information during data processing and which argparse
option flags to use when specifying them in python scripts. The `10x Genomics barcoded BAM format
<https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bam>`_ is followed
where that information is defined.

..  csv-table::
    :header: "SAM-tag", "Option flag", "Description"
    :widths: 20, 20, 40

    "BX", "-b, --barcode-tag", "SAM tag for storing the error corrected barcode."
    "MI", "-m, --molecule-tag", "SAM tag for storing molecule index specifying a identified molecule
    for each barcode."
    "MN", "-n, --number-tag", "SAM tag for storing molecule count for a particular barcode."
    "RX", "-s, --sequence-tag", "SAM tag for storing the uncorrected barcode sequence."
