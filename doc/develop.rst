Development
===========

Conda environment files
-----------------------

There are two types of files that describe Conda environments.

- The file ``environment.yml`` contains abstract dependencies such as ``pysam``,
  ``bowtie2`` or ``picard=2.10``. This file is managed manually and needs to be
  updated whenever there are new dependencies or when the required version for a
  dependency changes.

- The ``environment.linux.lock.yml`` and ``environment.osx.lock.yml`` files
  (lock files) contain a fully specified description of the entire environment,
  with locked-down versions.  These files are used to create the test
  environment.

Use the script ``misc/condalock.sh`` to update the lock files whenever you make
changes to ``environment.yml``.
