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

Then, you may need to edit the configuration file `workdir/blr.yaml`, in
particular to enter the path to your reference genome.

Finally, change into the `workdir` folder and run the pipeline:

    blr run


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

## Old version

To run the analysis described in [High throughput barcoding method for genome-scale phasing](https://www.nature.com/articles/s41598-019-54446-x),
look at the [stable branch](https://github.com/FrickTobias/BLR/tree/stable) for this git repository.

That version of BLR Analysis is also available at [OMICtools](https://omictools.com/blr-tool).
