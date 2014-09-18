import glob
import shutil
import os
import sys
import subprocess
import tarfile
import argparse
import vcf

import utils


class DestructWrapper(object):

    def __init__(self, install_directory):

        self.config_filename = os.path.join(install_directory, 'config.ini')
        self.destruct_script = os.path.join(os.path.dirname(__file__), os.path.pardir, 'destruct.py')


    def install(self):

        pass


    def run(self, temp_directory, bam_filenames, output_filename):

        bam_list_filename = os.path.join(temp_directory, 'bam_list.tsv')

        with open(bam_list_filename, 'w') as bam_list_file:
            for lib_id, bam_filename in bam_filenames.iteritems():
                bam_list_file.write(lib_id + '\t' + bam_filename + '\n')

        breakpoints_filename = output_filename
        breakreads_filename = os.path.join(temp_directory, 'breakreads.tsv')
        plots_tar_filename = os.path.join(temp_directory, 'plots.tar')
        destruct_tmp_directory = os.path.join(temp_directory, 'tmp')

        destruct_cmd = list()
        destruct_cmd += [sys.executable]
        destruct_cmd += [self.destruct_script]
        destruct_cmd += [bam_list_filename]
        destruct_cmd += [breakpoints_filename]
        destruct_cmd += [breakreads_filename]
        destruct_cmd += [plots_tar_filanem]
        destruct_cmd += ['--config', self.config_filename]
        destruct_cmd += ['--tmp', destruct_tmp_directory]
        destruct_cmd += ['--nocleanup', '--repopulate', '--maxjobs', 4, '--verbose']

        subprocess.check_call(destruct_cmd)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('install_directory', help='Destruct installation directory')
    args = parser.parse_args()

    delly = DellyWrapper(args.install_directory)

    delly.install()





