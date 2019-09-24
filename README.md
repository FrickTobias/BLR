[![Travis CI](https://api.travis-ci.org/FrickTobias/BLR.svg?branch=master)](https://travis-ci.org/FrickTobias/BLR/)

:exclamation:**NB! This is currently under heavy development.**:exclamation:

# Barcode-Linked Reads Analysis

## Usage

After installation, activate your Conda environment:

    conda activate blr

Choose a name for the analysis. It will be `workdir` in this example. Create
the analysis directory with this command:

    blr init --reads1=path/to/sample.R1.fastq.gz workdir

Note that BLR expects paired-end reads. However, only the path to the R1 file
needs to be provided. The R2 file will be found automatically.

Then, you may need to edit the configuration file `workdir/config.yaml`, in
particular to enter the path to your reference genome.

Finally, change into the `workdir` folder and run the pipeline:

    blr run reads.1.final.fastq.gz reads.2.final.fastq.gz


## Install and Setup

One-time installation:
- Clone the current repository

        git clone https://github.com/FrickTobias/BLR.git

- [Install miniconda](https://docs.conda.io/en/latest/miniconda.html)
- Enable the [bioconda channel](http://bioconda.github.io/)
- Create a new Conda environment with the dependencies:

      conda env create -n blr -f environment.yml

- Install blr into the environment in "editable install" mode:

      conda activate blr
      pip install -e .

This will install blr in such a way that you can still modify the source code
and get any changes immediately without re-installing.

## Pre-print version

To run the analysis described in [Efficient whole genome haplotyping and 
high-throughput single molecule phasing with barcode-linked reads](https://www.biorxiv.org/content/early/2018/06/26/356121)
look at the [stable branch](https://github.com/FrickTobias/BLR/tree/stable) for this git repository.

BLR Analysis is also available at [OMICtools](https://omictools.com/blr-tool).
