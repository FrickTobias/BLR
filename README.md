[![Travis CI](https://api.travis-ci.org/FrickTobias/BLR.svg)](https://travis-ci.org/FrickTobias/BLR/)

:exclamation:**NB! This is currently under heavy development.**:exclamation:

# Barcode-Linked Reads Analysis

## Usage
Create a working directory in which to run the analysis. 

    mkdir workdir

Create soft links for the gzipped fastq files to analyze in the `workdir` folder for read 1 and 2. Note that the  must have the names 
`reads.1.fastq.gz` and `reads.2.fastq.gz`

    ln -s $PWD/path/to/sample.R1.fastq.gz workdir/reads.1.fastq.gz
    ln -s $PWD/path/to/sample.R2.fastq.gz workdir/reads.2.fastq.gz

Make a copy of the `config.yaml` file and enter the path to your reference genome.
  
    cp config.yaml workdir/my_config.yaml
    
Activate your enviroment

    conda activate blr

Run the pipeline using snakemake and indicate your intended target. In this example the full pipeline is run.  

    snakemake --snakefile path/to/Snakefile --configfile workdir/my_config.yaml workdir/reads.1.final.fastq.gz workdir/reads.2.final.fastq.gz 

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
