import glob
import shutil
import os
import sys
import subprocess
import tarfile
import argparse
import vcf
import pandas as pd

import utils


class DellyWrapper(object):

    features = ['num_spanning', 'num_split']

    def __init__(self, install_directory):

        self.install_directory = os.path.abspath(install_directory)

        self.packages_directory = os.path.join(self.install_directory, 'packages')
        self.bin_directory = os.path.join(self.install_directory, 'bin')
        self.data_directory = os.path.join(self.install_directory, 'data')

        self.bamtools_root = os.path.join(self.packages_directory, 'bamtools')
        self.seqtk_root = os.path.join(self.packages_directory, 'seqtk')
        self.boost_root = self.install_directory

        self.boost_url = 'http://downloads.sourceforge.net/project/boost/boost/1.56.0/boost_1_56_0.tar.gz'
        self.boost_package = 'boost_1_56_0'
        self.boost_tgz = self.boost_package+'.tar.gz'

        self.delly_bin = os.path.join(self.bin_directory, 'delly')
        self.delly_excl_chrom = os.path.join(self.packages_directory, 'delly', 'human.hg19.excl.tsv')

        ensembl_version = '70'
        ensembl_genome_version = 'GRCh37'
        self.primary_assembly_url = 'ftp://ftp.ensembl.org/pub/release-'+ensembl_version+'/fasta/homo_sapiens/dna/Homo_sapiens.'+ensembl_genome_version+'.'+ensembl_version+'.dna.primary_assembly.fa.gz'
        self.chromosome_url = 'ftp://ftp.ensembl.org/pub/release-'+ensembl_version+'/fasta/homo_sapiens/dna/Homo_sapiens.'+ensembl_genome_version+'.'+ensembl_version+'.dna.chromosome.{0}.fa.gz'
        self.genome_fasta = os.path.join(self.data_directory, 'genome.fa')


    def install(self, **kwargs):

        Sentinal = utils.SentinalFactory(os.path.join(self.install_directory, 'sentinal_'), kwargs)

        with Sentinal('install_bamtools') as sentinal:

            if sentinal.unfinished:

                with utils.CurrentDirectory(self.packages_directory):

                    utils.rmtree('bamtools')
                    subprocess.check_call('git clone https://github.com/pezmaster31/bamtools.git', shell=True)

                    with utils.CurrentDirectory('bamtools'):

                        utils.makedirs('build')

                        with utils.CurrentDirectory('build'):

                            subprocess.check_call('cmake ..', shell=True)
                            subprocess.check_call('make', shell=True)


        with Sentinal('install_kseq') as sentinal:

            if sentinal.unfinished:

                with utils.CurrentDirectory(self.packages_directory):

                    utils.rmtree('seqtk')
                    subprocess.check_call('git clone https://github.com/lh3/seqtk.git', shell=True)


        with Sentinal('install_boost') as sentinal:

            if sentinal.unfinished:

                with utils.CurrentDirectory(self.packages_directory):

                    subprocess.check_call('wget ' + self.boost_url, shell=True)
                    subprocess.check_call('tar -xzvf ' + self.boost_tgz, shell=True)

                    with utils.CurrentDirectory(self.boost_package):

                        subprocess.check_call('./bootstrap.sh', shell=True)

                        clean_env = 'INCLUDE_PATH= CPLUS_INCLUDE_PATH='

                        boost_build_cmd = clean_env + ' ./b2 install'
                        boost_build_cmd += ' link=static'
                        boost_build_cmd += ' --prefix=' + self.install_directory
                        boost_build_cmd += ' --with-iostreams'
                        boost_build_cmd += ' --with-filesystem'
                        boost_build_cmd += ' --with-program_options'
                        boost_build_cmd += ' --with-date_time'

                        subprocess.check_call(boost_build_cmd, shell=True)


        with Sentinal('download_delly') as sentinal:

            if sentinal.unfinished:

                with utils.CurrentDirectory(self.packages_directory):

                    utils.rmtree('delly')
                    subprocess.check_call('git clone https://github.com/tobiasrausch/delly.git', shell=True)

                    with utils.CurrentDirectory('delly'):

                        with open('Makefile.tmp', 'w') as f:
                            subprocess.check_call('sed s/-O9/-g\ -O3/g Makefile', shell=True, stdout=f)
                        os.rename('Makefile.tmp', 'Makefile')


        with Sentinal('install_delly') as sentinal:

            if sentinal.unfinished:

                with utils.CurrentDirectory(os.path.join(self.packages_directory, 'delly')):

                    make_cmd = 'make -B'
                    make_cmd += ' BAMTOOLS_ROOT=' + self.bamtools_root
                    make_cmd += ' SEQTK_ROOT=' + self.seqtk_root
                    make_cmd += ' BOOST_ROOT=' + self.boost_root
                    make_cmd += ' src/delly'

                    subprocess.check_call(make_cmd, shell=True)

            utils.symlink(os.path.join(self.packages_directory, 'delly', 'src', 'delly'),
                link_directory=self.bin_directory)


        with Sentinal('download_genome') as sentinal:

            if sentinal.unfinished:

                with utils.CurrentDirectory(os.path.join(self.data_directory)):

                    if kwargs.get('chromosomes', None) is None:

                        utils.wget_file_gunzip(self.primary_assembly_url, self.genome_fasta)

                    else:

                        with open(self.genome_fasta, 'w') as genome_file:

                            for chromosome in kwargs['chromosomes']:

                                chromosome_filename = './chromosome_{0}.fa'.format(chromosome)

                                utils.wget_file_gunzip(self.chromosome_url.format(chromosome), chromosome_filename)

                                with open(chromosome_filename, 'r') as chromosome_file:
                                    shutil.copyfileobj(chromosome_file, genome_file)

                                os.remove(chromosome_filename)


    def run(self, bam_filenames, output_filename, temp_directory, control_id=None):

        utils.makedirs(temp_directory)

        bams = list()
        for lib_id, bam_filename in bam_filenames.iteritems():
            bams += [utils.symlink(bam_filename, link_name='{0}.bam'.format(lib_id), link_directory=temp_directory)]
            utils.symlink(bam_filename+'.bai', link_name='{0}.bam.bai'.format(lib_id), link_directory=temp_directory)

        table_filenames = list()

        for sv_type in ('DEL', 'DUP', 'INV', 'TRA'):

            vcf_filename = os.path.join(temp_directory, 'delly_{0}.vcf'.format(sv_type))
            table_filename = os.path.join(temp_directory, 'delly_{0}.tsv'.format(sv_type))

            self.run_sv_type(sv_type, bams, vcf_filename)

            if not os.path.exists(vcf_filename):
                continue

            self.convert_output(vcf_filename, table_filename, control_id=control_id)

            table_filenames.append(table_filename)

        self.merge_tables(table_filenames, output_filename)


    def run_sv_type(self, sv_type, bams, output_filename):

        delly_cmd = list()
        delly_cmd += [self.delly_bin]
        delly_cmd += ['-t', sv_type]
        delly_cmd += ['-x', self.delly_excl_chrom]
        delly_cmd += ['-o', output_filename]
        delly_cmd += ['-g', self.genome_fasta]
        delly_cmd += bams

        subprocess.check_call(delly_cmd)


    def convert_output(self, vcf_filename, table_filename, control_id=None):

        vcf_reader = vcf.Reader(filename=vcf_filename)

        breakpoint_table = list()
        counts_table = list()

        for row in vcf_reader:

            prediction_id = row.ID

            chrom_1 = row.CHROM
            chrom_2 = row.INFO['CHR2']

            strand_1, strand_2 = [('-', '+')[a == '3'] for a in row.INFO['CT'].split('to')]

            coord_1 = row.POS
            coord_2 = row.sv_end

            if 'LowQual' in row.FILTER:
                qual = 0
            else:
                qual = 1

            breakpoint_table.append((prediction_id, chrom_1, chrom_2, strand_1, strand_2, coord_1, coord_2, qual))

            for call in row.samples:

                library = call.sample

                num_spanning = call.data.DV
                num_split = call.data.RV

                counts_table.append((prediction_id, library, num_spanning, num_split))

        breakpoint_table = pd.DataFrame(breakpoint_table, columns=['prediction_id',
                                                                   'chromosome_1', 'chromosome_2',
                                                                   'strand_1', 'strand_2',
                                                                   'position_1', 'position_2',
                                                                   'qual'])

        counts_table = pd.DataFrame(counts_table, columns=['prediction_id', 'library', 'num_spanning', 'num_split'])

        spanning_read_counts = counts_table.groupby('prediction_id')['num_spanning'].sum().reset_index()
        split_read_counts = counts_table.groupby('prediction_id')['num_split'].sum().reset_index()

        total_read_counts = counts_table.set_index(['prediction_id', 'library']).sum(axis=1).unstack()
        total_read_counts.columns = [a + '_count' for a in total_read_counts.columns]
        total_read_counts = total_read_counts.reset_index()

        breakpoint_table = breakpoint_table.merge(spanning_read_counts, on='prediction_id')
        breakpoint_table = breakpoint_table.merge(split_read_counts, on='prediction_id')
        breakpoint_table = breakpoint_table.merge(total_read_counts, on='prediction_id')

        # Filter based on evidence in control dataset
        if control_id is not None:
             breakpoint_table = breakpoint_table[breakpoint_table['{0}_count'.format(control_id)] == 0]          

        breakpoint_table.to_csv(table_filename, sep='\t', index=False)


    def merge_tables(self, input_filenames, output_filename):

        output_table = list()

        for input_filename in input_filenames:
            output_table.append(pd.read_csv(input_filename, sep='\t'))

        output_table = pd.concat(output_table, ignore_index=True)

        output_table.to_csv(output_filename, sep='\t', index=False)


if __name__ == '__main__':

    cmdline.interface(DellyWrapper)





