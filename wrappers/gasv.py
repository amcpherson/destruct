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


class GASVWrapper(object):

    features = ['tumour_count', 'num_split']

    def __init__(self, install_directory):

        self.install_directory = os.path.abspath(install_directory)

        self.packages_directory = os.path.join(self.install_directory, 'packages')
        self.bin_directory = os.path.join(self.install_directory, 'bin')
        self.data_directory = os.path.join(self.install_directory, 'data')

        self.bamtogasv_jar = os.path.join(self.bin_directory, 'BAMToGASV.jar')
        self.gasv_jar = os.path.join(self.bin_directory, 'GASV.jar')


    def install(self, **kwargs):

        Sentinal = utils.SentinalFactory(os.path.join(self.install_directory, 'sentinal_'), kwargs)

        with Sentinal('download_gasv') as sentinal:

            if sentinal.unfinished:

                with utils.CurrentDirectory(self.packages_directory):

                    utils.rmtree('gasv')
                    subprocess.check_call('svn checkout http://gasv.googlecode.com/svn/trunk/ gasv', shell=True)

        with Sentinal('install_gasv') as sentinal:

            if sentinal.unfinished:

                with utils.CurrentDirectory(os.path.join(self.packages_directory, 'gasv')):

                    subprocess.check_call('bash install', shell=True)

                with utils.CurrentDirectory(self.bin_directory):

                    utils.symlink(os.path.join(self.packages_directory, 'gasv', 'bin', 'BAMToGASV.jar'))
                    utils.symlink(os.path.join(self.packages_directory, 'gasv', 'bin', 'GASV.jar'))


    def run(self, temp_directory, bam_filenames, output_filename):

        temp_directory = os.path.abspath(temp_directory)

        utils.makedirs(temp_directory)

        sample_ids = list()
        bam_linknames = list()

        for sample_id, bam_filename in bam_filenames.iteritems():

            bam_linkname = os.path.join(temp_directory, 'sample_{0}.bam'.format(sample_id))

            utils.remove(bam_linkname)
            os.symlink(bam_filename, bam_linkname)

            bam_linknames.append(bam_linkname)

        header_filename = os.path.join(temp_directory, 'merged.header')

        with open(header_filename, 'w'):
            for sample_id in sample_ids:
                header_file.write('@RG\tID:{0}\tSM:{0}\tLB:{0}\tPL:Illumina\n'.format(sample_id))

        merged_bam_filename = os.path.join(temp_directory, 'merged.bam')

        merge_cmd = ' '.join(['samtools', 'merge', '-rh', header_filename, '-'] + bam_linknames)

        with open(merged_bam_filename, 'wb') as merge_bam_file:
            subprocess.check_call(merge_cmd, shell=True, stdout=merge_bam_file)

        gasv_files_prefix = os.path.join(temp_directory, 'gasv_inputs')

        bamtogasv_cmd = 'java -Xms512m -Xmx2048m -jar {0} {1} -OUTPUT_PREFIX {2} -LIBRARY_SEPARATED sep'.format(
            self.bamtogasv_jar, merged_bam_filename, gasv_files_prefix)

        print bamtogasv_cmd
        subprocess.check_call(bamtogasv_cmd, shell=True)

        gasv_pr_filename = gasv_files_prefix + '.gasv.in'
        # gasv_info_filename = gasv_files_prefix + '.info'

        # gasv_info = pd.read_csv(gasv_info_filename, sep='\t').iloc[0]

        # batches.append((gasv_pr_filename, 'PR', gasv_info['Lmin'].astype(str), gasv_info['Lmax'].astype(str)))

        # batch_filename = os.path.join(temp_directory, 'batch.in')

        # with open(batch_filename, 'w') as batch_file:
        #     for info in batches:
        #         batch_file.write('\t'.join(info) + '\n')

        gasv_cmd = 'java -jar {0} --batch {1}'.format(
            self.gasv_jar, gasv_pr_filename)
        print gasv_cmd

        subprocess.check_call(gasv_cmd, shell=True)


    def convert_output(self, vcf_filename, bam_filenames, table_filename):

        pass


    def merge_tables(self, input_filenames, output_filename):

        pass


if __name__ == '__main__':

    cmdline.interface(GASVWrapper)

