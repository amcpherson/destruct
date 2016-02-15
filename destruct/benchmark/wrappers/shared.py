import os
import shutil

import utils

ensembl_version = '70'
ensembl_genome_version = 'GRCh37'

primary_assembly_url = 'ftp://ftp.ensembl.org/pub/release-'+ensembl_version+'/fasta/homo_sapiens/dna/Homo_sapiens.'+ensembl_genome_version+'.'+ensembl_version+'.dna.primary_assembly.fa.gz'

chromosome_url = 'ftp://ftp.ensembl.org/pub/release-'+ensembl_version+'/fasta/homo_sapiens/dna/Homo_sapiens.'+ensembl_genome_version+'.'+ensembl_version+'.dna.chromosome.{0}.fa.gz'


def download_genome(genome_fasta, chromosomes=None):

    if chromosomes is None:

        utils.wget_file_gunzip(primary_assembly_url, genome_fasta)

    else:

        with open(genome_fasta, 'w') as genome_file:

            for chromosome in chromosomes:

                chromosome_filename = genome_fasta + '.chromosome_{0}.fa'.format(chromosome)

                utils.wget_file_gunzip(chromosome_url.format(chromosome), chromosome_filename)

                with open(chromosome_filename, 'r') as chromosome_file:
                    shutil.copyfileobj(chromosome_file, genome_file)

                os.remove(chromosome_filename)

