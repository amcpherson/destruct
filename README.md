# Destruct

Destruct is a tool for joint prediction of rearrangement breakpoints from single or multiple tumour samples.

## Installation

### Installing from conda

The recommended method of installation for destruct is using `conda`.  First install [anaconda python](https://store.continuum.io/cshop/anaconda/) from the continuum website.  Then add my channel, and the bioconda channel, and install destruct as follows.

    conda config --add channels https://conda.anaconda.org/dranew
    conda config --add channels 'bioconda'
    conda install destruct

### Installing from source

#### Clone Source Code

To install the code, first clone from bitbucket.  A recursive clone is preferred to pull in all submodules.

    git clone --recursive git@bitbucket.org:dranew/destruct.git

#### Dependencies

To install from source you will need several dependencies.  A list of dependencies can be found in the `conda` `yaml` file in the repo at `conda/destruct/meta.yaml`.

#### Build executables and install

To build executables and install the destruct code as a python package run the following command in the destruct repo:

    python setup.py install

## Setup

Download and setup of the reference genome is automated.  Select a directory on your system that will contain the reference data, herein referred to as `$ref_data_dir`.  The `$ref_data_dir` directory will be used in many of the subsequent scripts when running destruct.

Download the reference data and build the required indexes:

    destruct create_ref_data $ref_data_dir

## Run

### Input Data

Destruct takes multiple bam files as input.  Bam files can be multiple samples from the same patient, or a cohort of patients.

### Running Destruct

Running destruct involves invoking a single command, `destruct run`.  The result of destruct is a set of tables in TSV format.  Suppose we would like to run destruct on tumour bam file `$tumour_bam` and normal bam file `$normal_bam`.  The following command will predict breakpoints jointly on these bam files:

    destruct run $ref_data_dir \
        $breakpoint_table $breakpoint_library_table \
        $breakpoint_read_table \
        --bam_files $tumour_bam $normal_bam \
        --lib_ids tumour normal \
        --tmpdir $tmp

where `$breakpoint_table`, `$breakpoint_library_table`, and `$breakpoint_read_table ` are output tables, and `$tmp` is a unique temporary directory.  If you need to stop and restart the script, using the same temporary directory will allow the scripts to restart where it left off.

For parallelism options see the section [Parallelism using pypeliner](#markdown-header-parallelism-using-pypeliner).

## Testing Your Installation

To test your install, a script is provided to generate a bam.  First install dwgsim.

    conda install dwgsim bwa

Clone the destruct repo to obtain additional scripts.  Change to the destruct repo directory.

    git clone git@bitbucket.org:dranew/destruct.git
    cd destruct

Setup a chromosome 20, 21 specific reference dataset.

    destruct create_ref_data \
        -c examples/chromosome_20_21_user_config.py \
        destruct_ref_data/

Generate a bam file using the dwgsim read simulator.

    python destruct/benchmark/generate_bam.py \
        examples/breakpoint_simulation_config.yaml \
        destruct_ref_data/ \
        test_raw_data \
        test.bam \
        test.fasta \
        test_breakpoints.tsv \
        --submit local

Run destruct on the simulated bam file.

    destruct run \
        --config examples/chromosome_20_user_config.py \
        destruct_ref_data/ \
        breaks.tsv break_libs.tsv break_reads.tsv \
        --bam_files test.bam \
        --lib_ids test_sample \
        --raw_data_dir destruct_raw_data/ \
        --submit local

You can then compare the output breaks.tsv to test_breakpoints.tsv.

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
* `num_unique_reads`: Total number of discordant reads, potential PCR duplicates removed
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
* `template_length_min`: minimum of `template_length_1` and `template_length_2`

#### Breakpoint Library Table

The breakpoint library table contains information per breakpoint per input dataset.  The file is tab separated with the first line as the header and contains the following columns.

* `prediction_id`: Unique identifier of the breakpoint prediction
* `library_id`: ID of the dataset as given on the command line or in the input dataset table
* `num_reads`: Number of reads for this dataset
* `num_unique_reads`: Number of reads for this dataset, potential PCR duplicates removed

#### Breakpoint Read Table

The breakpoint library table contains information per breakpoint per discordant read.  The file is tab separated with the first line as the header and contains the following columns.

* `prediction_id`: Unique identifier of the breakpoint prediction
* `library_id`: ID of the dataset as given on the command line or in the input dataset table
* `fragment_id`: ID of the discordant read
* `read_end`: End of the paired end read
* `seq`: Read sequence
* `qual`: Read quality
* `comment`: Read comment
* `filtered`: The read was filtered prior to finalizing the prediction

## Parallelism Using Pypeliner

Destruct uses the pypeliner python library for parallelism.  Several of the scripts described above will complete more quickly on a multi-core machine or on a cluster.

To run a script in multicore mode, using a maximum of 4 cpus, add the following command line option:

    --maxjobs 4

To run a script on a cluster with qsub/qstat, add the following command line option:

    --submit asyncqsub

Often a call to qsub requires specific command line parameters to request the correct queue, and importantly to request the correct amount of memory.  To allow correct calls to qsub, use the `--nativespec` command line option, and use the placeholder `{mem}` which will be replaced by the amount of memory (in gigabytes) required for each job launched with qsub.  For example, to use qsub, and request queue `all.q` and set the `mem_free` to the required memory, add the following command line options:

    --submit asyncqsub --nativespec "-q all.q -l mem_free={mem}G"

# Build

## Docker builds

To build a destruct docker image, for instance version v0.4.13, run the following docker command:

    docker build --build-arg app_version=0.4.17 -t amcpherson/destruct:0.4.17 .
    docker push amcpherson/destruct:0.4.17

# Benchmarking

## Setup

The following requirements should be installed with conda:

    conda create -n lumpy lumpy-sv
    conda create -n delly -c dranew delly==0.7.3
    conda create -n sambamba sambamba
    conda create -n bcftools bcftools
    conda create -n vcftools -c bioconda perl-vcftools-vcf
    conda create -n htslib htslib

Add the requirements to your path with:

    PATH=$PATH:/opt/conda/envs/lumpy/bin/
    PATH=$PATH:/opt/conda/envs/delly/bin/
    LD_LIBRARY_PATH=/opt/conda/envs/delly/lib/
    PATH=$PATH:/opt/conda/envs/sambamba/bin/
    PATH=$PATH:/opt/conda/envs/bcftools/bin/
    PATH=$PATH:/opt/conda/envs/vcftools/bin/
    PATH=$PATH:/opt/conda/envs/htslib/bin/

You may also need to link sambamba into the lumpy install:

    ln -s /opt/conda/envs/sambamba/bin/sambamba /opt/conda/envs/lumpy/bin/sambamba

