![Travis CI badge](https://api.travis-ci.org/FrickTobias/BLR.svg)

# Barcode-Linked Reads Analysis

The blr analysis is currently under development.

## Usage

    bash BLR_automation [options] <read1> <read2> <output>

## Setup

One-time installation:

- [Install miniconda](https://docs.conda.io/en/latest/miniconda.html)
- Enable the [bioconda channel](http://bioconda.github.io/)
- Create a new Conda environment with the dependencies:

      conda env create -n blr -f environment.yml

- Install blr into the environment in "editable install" mode:

      conda activate blr
      pip install -e .

This will install blr in such a way that you can still modify the source code
and get any changes immediately without re-installing.

To run the program, run `blr`.

Subsequently, you will only need to activate the environment with

    conda activate blr

## Pre-print version

To run the analysis described in [Efficient whole genome haplotyping and 
high-throughput single molecule phasing with barcode-linked reads](https://www.biorxiv.org/content/early/2018/06/26/356121)
look at the [stable branch](https://github.com/FrickTobias/BLR/tree/stable) for this git repository.

BLR Analysis is also available at [OMICtools](https://omictools.com/blr-tool).
