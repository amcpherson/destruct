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


class LumpySVWrapper(object):

    features = ['tumour_count', 'num_split']

    def __init__(self, install_directory):

        self.install_directory = os.path.abspath(install_directory)

        self.packages_directory = os.path.join(self.install_directory, 'packages')
        self.bin_directory = os.path.join(self.install_directory, 'bin')
        self.data_directory = os.path.join(self.install_directory, 'data')

        self.lumpysv_bin = os.path.join(self.bin_directory, 'lumpysv')


    def install(self, **kwargs):

        Sentinal = utils.SentinalFactory(os.path.join(self.install_directory, 'sentinal_'), kwargs)

        with Sentinal('download_lumpysv') as sentinal:

            if sentinal.unfinished:

                with utils.CurrentDirectory(self.packages_directory):

                    utils.rmtree('lumpy-sv')
                    subprocess.check_call('git clone git://github.com/arq5x/lumpy-sv.git', shell=True)

        with Sentinal('install_lumpysv') as sentinal:

            if sentinal.unfinished:

                with utils.CurrentDirectory(os.path.join(self.packages_directory, 'lumpy-sv')):

                    subprocess.check_call('make', shell=True)

                    utils.makedirs(self.bin_directory)
                    utils.symlink(os.path.join('bin', 'lumpy'), link_directory=self.bin_directory)


    def run(self, temp_directory, bam_filenames, output_filename):

        pass


    def convert_output(self, vcf_filename, bam_filenames, table_filename):

        pass


    def merge_tables(self, input_filenames, output_filename):

        pass


if __name__ == '__main__':

    cmdline.interface(LumpySVWrapper)



