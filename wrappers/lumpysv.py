import glob
import shutil
import os
import sys
import subprocess
import tarfile
import argparse
import vcf
import numpy as np
import pandas as pd

import shared
import utils
import cmdline
import pypeliner.commandline


class LumpySVWrapper(object):

    features = ['evidence_set_score', 'read_count']

    def __init__(self, install_directory):

        self.wrappers_bin_directory = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'bin')

        self.install_directory = os.path.abspath(install_directory)

        self.packages_directory = os.path.join(self.install_directory, 'packages')
        self.bin_directory = os.path.join(self.install_directory, 'bin')
        self.data_directory = os.path.join(self.install_directory, 'data')

        self.lumpy_bin = os.path.join(self.bin_directory, 'lumpy')
        self.scripts_directory = os.path.join(self.packages_directory, 'lumpy-sv', 'scripts')

        self.genome_fasta = os.path.join(self.data_directory, 'genome.fa')


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

        with Sentinal('install_yaya') as sentinal:

            if sentinal.unfinished:

                with utils.CurrentDirectory(self.packages_directory):

                    utils.rmtree('yaha')
                    subprocess.check_call('git clone https://github.com/GregoryFaust/yaha.git', shell=True)

        with Sentinal('download_genome') as sentinal:

            if sentinal.unfinished:

                utils.makedirs(self.data_directory)
                shared.download_genome(self.genome_fasta, kwargs.get('chromosomes', None))

        with Sentinal('index_genome') as sentinal:

            if sentinal.unfinished:

                subprocess.check_call('bwa index '+self.genome_fasta, shell=True)


    def prepare(self, bam, prefix, pe_id, sr_id):

        unmapped_fastq = prefix + 'um.fastq'

        pypeliner.commandline.execute(
            *'samtools view {paired_bam} \
            | {scripts_directory}/split_unmapped_to_fasta.pl -b 20 \
            > {unmapped_fastq}' \
                .format(
                    paired_bam=bam, 
                    scripts_directory=self.scripts_directory,
                    unmapped_fastq=unmapped_fastq) \
                .split())

        split_unsort_bam = prefix + 'split.unsort.bam'

        pypeliner.commandline.execute(
            *'bwa bwasw -H -t 20 {genome} {unmapped_fastq} \
            | samtools view -Sb - > {split_unsort_bam}' \
                .format(
                    genome=self.genome_fasta,
                    unmapped_fastq=unmapped_fastq,
                    split_unsort_bam=split_unsort_bam) \
                .split())

        split_prefix = prefix + 'split'
        split_bam = split_prefix + '.bam'

        pypeliner.commandline.execute(
            *'samtools sort {split_unsort_bam} {split_prefix}' \
                .format(
                    split_unsort_bam=split_unsort_bam,
                    split_prefix=split_prefix) \
                .split())

        paired_end_histo = prefix + 'pe.histo'
        paired_end_dist = prefix + 'pe.dist'

        pypeliner.commandline.execute(
            *'samtools view {paired_bam} \
            | python {wrappers_bin_directory}/skip_lines.py 100000 \
            | {scripts_directory}/pairend_distro.py \
            -r 150 \
            -X 4 \
            -N 10000 \
            -o {paired_end_histo} \
            > {paired_end_dist}' \
                .format(
                    wrappers_bin_directory=self.wrappers_bin_directory,
                    scripts_directory=self.scripts_directory,
                    paired_bam=bam,
                    paired_end_histo=paired_end_histo,
                    paired_end_dist=paired_end_dist) \
                .split())

        read_length = prefix + 'pe.read_length'

        pypeliner.commandline.execute(
            *'samtools view {paired_bam} \
            | python {wrappers_bin_directory}/sam_max_read_length.py \
            > {read_length}' \
                .format(
                    wrappers_bin_directory=self.wrappers_bin_directory,
                    paired_bam=bam,
                    read_length=read_length
                    )
                .split()
            )

        pe_params = dict()

        with open(paired_end_dist, 'r') as f:
            for line in f:
                for entry in line.split():
                    key, value = entry.split(':')
                    pe_params[key] = float(value)

        with open(read_length, 'r') as f:
            pe_params['read_length'] = int(f.read())

        pe_params['bam_file'] = bam
        pe_params['histo_file'] = paired_end_histo
        pe_params['min_non_overlap'] = pe_params['read_length']
        pe_params['discordant_z'] = 4
        pe_params['back_distance'] = 20
        pe_params['weight'] = 1
        pe_params['min_mapping_threshold'] = 20
        pe_params['id'] = pe_id

        sr_params = dict()

        sr_params['bam_file'] = split_bam
        sr_params['back_distance'] = 20
        sr_params['weight'] = 1
        sr_params['min_mapping_threshold'] = 20
        sr_params['id'] = sr_id

        args = ''
        args += ' -pe ' + ','.join([':'.join((key, str(value))) for key, value in pe_params.iteritems()])
        args += ' -sr ' + ','.join([':'.join((key, str(value))) for key, value in sr_params.iteritems()])

        return args


    def run(self, tumour_bam, normal_bam, output_filename, temp_directory):

        utils.makedirs(temp_directory)

        sample_id_map = {
            1:'tumour_pe',
            2:'tumour_sr',
            3:'normal_pe',
            4:'normal_sr',
        }

        tumour_args = self.prepare(tumour_bam, os.path.join(temp_directory, 'tumour.'), 1, 2)

        normal_args = self.prepare(normal_bam, os.path.join(temp_directory, 'normal.'), 3, 4)

        results_bedpe = os.path.join(temp_directory, 'results.bedpe')

        pypeliner.commandline.execute(
            *'{lumpy_bin} \
            -t {lumpy_tmp} \
            -mw 2 \
            -tt 0.0 \
            {tumour_args} \
            {normal_args} \
            > {results_bedpe}' \
                .format(
                    lumpy_bin=self.lumpy_bin,
                    lumpy_tmp=os.path.join(temp_directory, 'lumpy.'),
                    tumour_args=tumour_args,
                    normal_args=normal_args,
                    results_bedpe=results_bedpe
                    )
                .split()
            )

        bedpe_columns = [
            'chromosome_1',
            'interval_1_start',
            'interval_1_end',
            'chromosome_2',
            'interval_2_start',
            'interval_2_end',
            'prediction_id',
            'evidence_set_score',
            'strand_1',
            'strand_2',
            'type',
            'id_of_samples',
            'strand_configurations',
            'breakpoint_max_prob',
            'segment_95_prob',
        ]

        results = pd.read_csv(results_bedpe, sep='\t', names=bedpe_columns,
            converters={'chromosome_1':str, 'chromosome_2':str})

        column_prefixes = [
            ('type', 'TYPE:'),
            ('id_of_samples', 'IDS:'),
            ('strand_configurations', 'STRANDS:'),
            ('breakpoint_max_prob', 'MAX:'),
            ('segment_95_prob', '95:'),
        ]

        for column, prefix in column_prefixes:
            results[column] = results[column].str.replace(prefix, '')

        results['position_1'] = None
        results['position_2'] = None

        results[['position_1', 'position_2']] = results['breakpoint_max_prob'].str.split(';').apply(pd.Series, 1)

        results['position_1'] = results['position_1'].str.split(':').apply(lambda a: a[1]).astype(int)
        results['position_2'] = results['position_2'].str.split(':').apply(lambda a: a[1]).astype(int)

        read_counts = list()

        for prediction_id, id_of_samples in results[['prediction_id', 'id_of_samples']].values:
            for sample_id, read_count in (a.split(',') for a in id_of_samples.split(';')):
                read_counts.append((prediction_id, int(sample_id), int(read_count)))

        read_counts = pd.DataFrame(read_counts, columns=('prediction_id', 'sample_id', 'read_count'))

        read_counts = read_counts.set_index(['prediction_id', 'sample_id'])['read_count']\
                                 .unstack()\
                                 .fillna(0)\
                                 .rename(columns=sample_id_map)

        # Filter germline
        read_counts = read_counts[(read_counts['normal_pe'] == 0) & (read_counts['normal_sr'] == 0)]

        read_counts = read_counts.rename(columns={'tumour_pe':'pe_read_count',
                                                  'tumour_sr':'sr_read_count'})

        read_counts['read_count'] = read_counts['pe_read_count'] + read_counts['sr_read_count']

        read_counts['normal_count'] = 0

        read_counts = read_counts.drop(['normal_pe', 'normal_sr'], axis=1)

        results = results.merge(read_counts, left_on='prediction_id', right_index=True)

        # Replace inf with large number
        results.loc[np.isposinf(results['evidence_set_score']), 'evidence_set_score'] = 1e6

        # Lower scores are better
        results['evidence_set_score'] = -results['evidence_set_score']

        results.to_csv(output_filename, sep='\t', index=False)


if __name__ == '__main__':

    cmdline.interface(LumpySVWrapper)



