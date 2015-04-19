# Destruct

Destruct is a tool for joint prediction of rearrangement breakpoints from single or multiple tumour samples.

## Prerequisites

### Python

Destruct requires python and the numpy/scipy stack.  The recommended way to install python (also the easiest), is to download and install the latest (Anaconda Python](https://store.continuum.io/cshop/anaconda/) from the Continuum website.

#### Python libraries

If you do no use anaconda, you will be required to install the following python libraries.

* [numpy/scipy](http://www.numpy.org)
* [pandas](http://pandas.pydata.org)
* [matplotlib](http://matplotlib.org)

### Scons

Building the source requires scons, verson 2.3.4 can be installed as follows:

    wget http://prdownloads.sourceforge.net/scons/scons-2.3.4.tar.gz
    tar -xzvf scons-2.3.4.tar.gz
    cd scons-2.3.4
    python setup.py install

### Samtools

[Samtools](http://www.htslib.org) is required and should be on the path.

### Bowtie

Destruct uses [Bowtie](http://bowtie-bio.sourceforge.net/index.html) for one of the alignment steps.  The `bowtie` and `bowtie-build` executables should be on the path.  Note that destruct does not use bowtie 2.

## Install

### Clone Source Code

To install the code, first clone from bitbucket.  A recursive clone is preferred to pull in all submodules.

    git clone --recursive git@bitbucket.org:dranew/destruct.git

The following steps will assume you are in the `destruct` directory.

    cd destruct

### Build Executables

Destruct requires compilation of a number of executables using scons.

    cd src
    scons install

### Install Python Libraries

There are two options for installing the python libraries.

#### Option 1:

A temporary solution is to modify the python path to point to the location of the source code on your system.  If you use bash, the following command will correctly modify your environment.

    source pythonpath.sh

#### Option 2:

A more permanent solution is to install the libraries into your python site packages.  Note that if you are using python installed on your system, you may need admin privileges.

To install pypeliner, a pipeline management system:

    cd pypeliner
    python setup.py install

## Setup

Download and setup of the reference genome is automated.  Select a directory on your system that will contain the reference data, herein referred to as `$ref_data_dir`.  The `$ref_data_dir` directory will be used in many of the subsequent scripts when running destruct.

Download the reference data and build the required indexes:

    python create_ref_data.py $ref_data_dir

## Run

### Input Data

Destruct takes multiple bam files as input.  Bam files can be multiple samples from the same patient, or a cohort of patients.

### Running Destruct

Running destruct involves invoking a single script, `destruct.py`.  The result of destruct is several tables and a tar of plots.  Suppose we would like to run destruct on tumour bam file `$tumour_bam` and normal bam file `$normal_bam`.  The following command will predict breakpoints jointly on these bam files:

    python destruct.py $ref_data_dir \
        $breakpoint_table $breakpoint_library_table $plots_tar \
        --breakpoint_read_table $breakpoint_read_table \
        --bam_files $tumour_bam $normal_bam \
        --lib_ids tumour normal \
        --tmpdir $tmp

where `$breakpoint_table`, `$breakpoint_library_table`, `$plots_tar`, and `$breakpoint_read_table ` are output tables and plots, and `$tmp` is a unique temporary directory.  If you need to stop and restart the script, using the same temporary directory will allow the scripts to restart where it left off.

For parallelism options see the section [Parallelism using pypeliner](#markdown-header-parallelism-using-pypeliner).

### Output File Formats

#### Breakpoint Table

The breakpoint library table contains information per breakpoint.  The file is tab separated with the first line as the header and contains the following columns.

* `prediction_id`: Unique identifier of the breakpoint prediction
* `chromosome_1`: Chromosome for breakend 1
* `strand_1`: Strand for breakend 1
* `position_1`: Position of breakend 1
* `gene_id1`: Ensembl gene id for gene at or near breakend 1
* `gene_name1`: Name of the gene at or near breakend 1
* `gene_location1`: Location of the gene with respect to the breakpoint for breakend 1
* `chromosome_2`: Chromosome for breakend 2
* `strand_2`: Strand for breakend 2
* `position_2`: Position of breakend 2
* `gene_id2`: Ensembl gene id for gene at or near breakend 2
* `gene_name2`: Name of the gene at or near breakend 2
* `gene_location2`: Location of the gene with respect to the breakpoint for breakend 2
* `type`: Breakpoint orientation type deletion: +-, inversion: ++ or --, duplication -+, translocation: 2 different chromosomes
* `num_reads`: Total number of discordant reads
* `num_split`: Total number of discordant reads split by the breakpoint
* `num_inserted`: Number of untemplated nucleotides inserted at the breakpoint
* `homology`: Sequence homology at the breakpoint
* `dgv_ids`: Database of genomic variants annotation for germline variants
* `sequence`: Sequence as predicted by discordant reads and possibly split reads
* `inserted`: Nucleotides inserted at the breakpoint
* `mate_score`: average score of mate reads aligning as if concordant
* `template_length_1`: length of region to which discordant reads align at breakend 1
* `template_length_2`: length of region to which discordant reads align at breakend 2
* `log_cdf`: mean cdf of discordant read alignment likelihood
* `log_likelihood`: mean likelihood of discordant read alignments
* `template_length_min`: minimum of `template_length_1` and `t`emplate_length_2`

#### Breakpoint Library Table

The breakpoint library table contains information per breakpoint per input dataset.  The file is tab separated with the first line as the header and contains the following columns.

* `prediction_id`: Unique identifier of the breakpoint prediction
* `library_id`: ID of the dataset as given on the command line or in the input dataset table
* `num_reads`: Number of reads for this dataset.

#### Plots Tar

The plots tar file contains pdfs of the score likelihoods per alignment length and fragment length distribution.

#### Breakpoint Read Table

The breakpoint library table contains information per breakpoint per discordant read.  The file is tab separated with the first line as the header and contains the following columns.

* `prediction_id`: Unique identifier of the breakpoint prediction
* `library_id`: ID of the dataset as given on the command line or in the input dataset table
* `fragment_id`: ID of the discordant read
* `read_end`: End of the paired end read 
* `seq`: Read sequence
* `qual`: Read quality
* `comment`: Read comment

## Parallelism Using Pypeliner

Destruct uses the pypeliner python library for parallelism.  Several of the scripts described above will complete more quickly on a multi-core machine or on a cluster.

To run a script in multicore mode, using a maximum of 4 cpus, add the following command line option:

    --maxjobs 4

To run a script on a cluster with qsub/qstat, add the following command line option:

    --submit asyncqsub

Often a call to qsub requires specific command line parameters to request the correct queue, and importantly to request the correct amount of memory.  To allow correct calls to qsub, use the `--nativespec` command line option, and use the placeholder `{mem}` which will be replaced by the amount of memory (in gigabytes) required for each job launched with qsub.  For example, to use qsub, and request queue `all.q` and set the `mem_free` to the required memory, add the following command line options:

    --submit asyncqsub --nativespec "-q all.q -l mem_free={mem}G"

momac01:destruct2 amcphers$