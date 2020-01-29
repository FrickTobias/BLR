[![Travis CI](https://api.travis-ci.org/FrickTobias/BLR.svg?branch=master)](https://travis-ci.org/FrickTobias/BLR/)

:exclamation:**NB! This is currently under heavy development.**:exclamation:

# Barcode-Linked Reads Analysis

Content:

- [Usage](#Usage)
- [One-time installation](#One-time-installation)
- [Old version](#Old-version)

## Usage

### 1. Setup an analysis folder

Activate your Conda environment.

    conda activate blr

Choose a name for the analysis. It will be `output_folder` in this example. Create
the analysis directory.

    blr init --reads1=path/to/sample.R1.fastq.gz path/to/output_folder

Note that BLR expects paired-end reads. However, only the path to the R1 file
needs to be provided. The R2 file will be found automatically.

To use the other blr commands, make sure you working directory is your 
newly created analysis folder.

    cd path/to/output_folder

Then, you may need to edit the configuration file `blr.yaml`, in
particular to enter the path to your reference genome. 

    blr set --genome_reference=path/to/GRCh38.fasta

To see what other configurations can be altered, read the documentation in 
the `blr.yaml` file.

### 2. Running an analysis

Change working directory to your analysis folder

    cd path/to/output_folder

The pipeline automatically runs all steps.

    blr run

For more options, see the documentation.

    blr run -h

## One-time installation

### 1. Prerequisite: Conda

- [Install miniconda](https://docs.conda.io/en/latest/miniconda.html)
- Enable the [bioconda channel](http://bioconda.github.io/)

You could also try copy-pasting the following to your terminal. This will download miniconda, 
install it to you `$HOME` folder and enable the bioconda channel.

    if [[ $OSTYPE = "linux-gnu" ]]; then 
        wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    elif [[ $OSTYPE = "darwin"* ]]; then 
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh
    fi
    bash miniconda.sh -b -p $HOME/miniconda
    source $HOME/miniconda/etc/profile.d/conda.sh
    conda config --add channels bioconda

### 2. Install

Clone the git repository.

      git clone https://github.com/FrickTobias/BLR.git

Create & activate a new Conda environment, in which all dependencies will be 
installed.

      conda env create -n blr -f environment.yml
      conda activate blr

Install blr into the environment in "editable install" mode.

      pip install -e .

This will install blr in such a way that you can still modify the source code
and get any changes immediately without re-installing.

#### 2.1 Linux users (not macOS)

To enable DeepVariant, install it separately to your environment.

    conda activate blr
    conda install deepvariant

This will enable the `variant_caller: deepvariant` option in the analysis config file.    

### 3. Updating

Change working directory to your blr git folder and update.

    cd path/to/BLR
    git pull

## Old version

To run the analysis described in [High throughput barcoding method for genome-scale phasing](https://www.nature.com/articles/s41598-019-54446-x),
look at the [stable branch](https://github.com/FrickTobias/BLR/tree/stable) for this git repository.

That version of BLR Analysis is also available at [OMICtools](https://omictools.com/blr-tool).