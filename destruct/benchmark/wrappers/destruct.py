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
import pypeliner
import cmdline


class DestructWrapper(object):

    features = ['num_reads', 'num_split', 'template_length_min', 'log_likelihood', 'log_cdf', 'mate_score']

    def __init__(self, install_directory):

        self.install_directory = os.path.abspath(install_directory)

        self.ref_data_directory = os.path.join(self.install_directory, 'data')
        self.user_config_filename = os.path.join(self.install_directory, 'user_config.py')

        self.createref_script = os.path.join(os.path.dirname(__file__), os.path.pardir, os.path.pardir, 'create_ref_data.py')
        self.destruct_script = os.path.join(os.path.dirname(__file__), os.path.pardir, os.path.pardir, 'run_destruct.py')


    def install(self, **kwargs):

        Sentinal = utils.SentinalFactory(os.path.join(self.install_directory, 'sentinal_'), kwargs)

        with Sentinal('createref') as sentinal:

            if sentinal.unfinished:

                utils.makedirs(self.install_directory)
                
                with open(self.user_config_filename, 'w') as user_config_file:
                    if kwargs.get('chromosomes', None) is not None:
                        chromosomes = kwargs['chromosomes']
                        ensembl_assemblies = ['chromosome.'+a for a in chromosomes]
                        user_config_file.write('chromosomes = '+repr(chromosomes)+'\n')
                        user_config_file.write('ensembl_assemblies = '+repr(ensembl_assemblies)+'\n')

                createref_cmd = [sys.executable]
                createref_cmd += [self.createref_script]
                createref_cmd += [self.ref_data_directory]
                createref_cmd += ['-c', self.user_config_filename]

                subprocess.check_call(createref_cmd)


    def run(self, bam_filenames, output_filename, temp_directory, control_id=None):

        utils.makedirs(temp_directory)

        lib_ids = bam_filenames.keys()
        bams = bam_filenames.values()

        breakpoint_table_filename = os.path.join(temp_directory, 'breakpoint.tsv')
        breakpoint_library_table_filename = os.path.join(temp_directory, 'breakpoint_library.tsv')
        destruct_tmp_directory = os.path.join(temp_directory, 'tmp')

        destruct_cmd = list()
        destruct_cmd += [sys.executable]
        destruct_cmd += [self.destruct_script]
        destruct_cmd += [self.ref_data_directory]
        destruct_cmd += [breakpoint_table_filename]
        destruct_cmd += [breakpoint_library_table_filename]

        destruct_cmd += ['--lib_ids']
        destruct_cmd += lib_ids

        destruct_cmd += ['--bam_files']
        destruct_cmd += bams

        destruct_cmd += ['--config', self.user_config_filename]
        destruct_cmd += ['--tmp', destruct_tmp_directory]
        destruct_cmd += ['--nocleanup', '--repopulate', '--maxjobs', '4', '--loglevel', 'DEBUG']

        pypeliner.commandline.execute(*destruct_cmd)

        breakpoint_table = pd.read_csv(breakpoint_table_filename, sep='\t')
        breakpoint_library_table = pd.read_csv(breakpoint_library_table_filename, sep='\t')

        breakpoint_counts = breakpoint_library_table.set_index(['prediction_id', 'library'])[['num_reads']]\
                                                    .unstack()\
                                                    .fillna(0)\
                                                    .astype(int)
        breakpoint_counts.columns = [a[1] + '_count' for a in breakpoint_counts.columns]

        breakpoint_table = breakpoint_table.merge(breakpoint_counts,
                                                  left_on='prediction_id',
                                                  right_index=True)
        
        # Filter based on evidence in control dataset
        if control_id is not None:
             breakpoint_table = breakpoint_table[breakpoint_table['{0}_count'.format(control_id)] == 0]          

        breakpoint_table.to_csv(output_filename, sep='\t', index=False)


if __name__ == '__main__':

    cmdline.interface(DestructWrapper)





