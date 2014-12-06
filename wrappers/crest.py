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
import cmdline


class CRESTWrapper(object):

    features = ['tumour_count', 'num_split']

    def __init__(self, install_directory):

        self.install_directory = install_directory

        self.packages_directory = os.path.join(self.install_directory, 'packages')
        self.data_directory = os.path.join(self.install_directory, 'data')
        self.crest_directory = os.path.join(self.packages_directory, 'crest')

        self.extractsclip_script = os.path.join(self.crest_directory, 'extractSClip.pl')
        self.countdiff_script = os.path.join(self.crest_directory, 'countDiff.pl')
        self.crest_script = os.path.join(self.crest_directory, 'CREST.pl')

	self.genome_fasta = os.path.join(self.data_directory, 


    def install(self, **kwargs):

        Sentinal = utils.SentinalFactory(os.path.join(self.install_directory, 'sentinal_'), kwargs)

        with Sentinal('download_crest') as sentinal:

            if sentinal.unfinished:

                with utils.CurrentDirectory(self.packages_directory):

                    utils.remove('CREST.tgz')
                    subprocess.check_call('wget ftp://ftp.stjude.org/pub/software/CREST/CREST.tgz', shell=True)

                    utils.makedirs('crest')
                    subprocess.check_call('tar -C crest -xzvf CREST.tgz', shell=True)


    def run(self, tumour_bam, normal_bam, output_filename, temp_directory):

        tumour_cover_filename = os.path.join(temp_directory, 'tumour.bam.cover')
        normal_cover_filename = os.path.join(temp_directory, 'normal.bam.cover')

        for cover_filename, bam in zip((tumour_cover_filename, normal_cover_filename),
                                       (tumour_bam, normal_bam)):

            cmd = ['perl', self.extractsclip_script]
            cmd += ['-i', tumour_bam]
            cmd += ['--ref_genome', self.genome_fasta]

            with open(cover_filename, 'w') as cover_file:
                subprocess.check_call(cmd, stdout=cover_file)

        cmd = ['perl', self.countdiff_script]
        cmd += ['-d', tumour_cover_filename]
        cmd += ['-g', normal_cover_filename]

        soft_clip_dist_filename = os.path.join(temp_directory, 'soft_clip.dist.txt')

        with open(soft_clip_dist_filename, 'w') as soft_clip_dist_file:
            subprocess.check_call(cmd, stdout=soft_clip_dist_file)

        somatic_cover_filename = tumour_cover_filename + '.somatic.cover'

        cmd = ['perl', self.crest_script]
        cmd += ['-f', somatic_cover_filename]
        cmd += ['-d', tumour_bam]
        cmd += ['-g', normal_bam]
        cmd += ['--ref_genome', self.genome_fasta]
        cmd += ['-t', self.genome_2bit]
        cmd += ['--blatport', self.blatport]


    def convert_output(self, vcf_filename, bam_filenames, table_filename):

        pass


    def merge_tables(self, input_filenames, output_filename):

        pass


if __name__ == '__main__':

    cmdline.interface(CRESTWrapper)



