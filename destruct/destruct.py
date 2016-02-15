import csv
import sys
import logging
import os
import ConfigParser
import re
import itertools
import subprocess
import argparse
import string
import tarfile
import collections
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import pygenes
import pypeliner
import pypeliner.workflow
import pypeliner.managed as mgd

import score_stats
import utils.plots
import utils.misc
import utils.seq
import predict_breaks


destruct_directory = os.environ.get('DESTRUCT_PACKAGE_DIRECTORY', None)
if destruct_directory is None:
    raise Exception('please set the $DESTRUCT_PACKAGE_DIRECTORY environment variable to the root of the destruct package')

data_directory = os.path.join(destruct_directory, 'data')
bin_directory = os.path.join(destruct_directory, 'bin')
default_config_filename = os.path.join(destruct_directory, 'defaultconfig.py')


__version__ = '0.1.0'

if __name__ == '__main__':

    import destruct

    argparser = argparse.ArgumentParser()

    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('--version', action='version', version=__version__)

    argparser.add_argument('ref_data_dir',
                           help='Reference dataset directory')

    argparser.add_argument('breakpoint_table',
                           help='Output table of breakpoint information in TSV format')

    argparser.add_argument('breakpoint_library_table',
                           help='Output table of library specific breakpoint information in TSV format')

    argparser.add_argument('plots_tar',
                           help='Output diagnostic plots tar filename')

    argparser.add_argument('--breakpoint_read_table', required=False,
                           help='Output table of breakpoint read information in TSV format')

    argparser.add_argument('--libs_table', required=False,
                           help='Input libraries list table filename')

    argparser.add_argument('--bam_files', nargs='+', required=False,
                           help='Input bam filenames')

    argparser.add_argument('--lib_ids', nargs='+', required=False,
                           help='Input ids for respective bam filenames')

    argparser.add_argument('--config', required=False,
                           help='Configuration Filename')

    args = vars(argparser.parse_args())

    if not ((args['libs_table'] is not None) or (args['bam_files'] is not None and args['lib_ids'] is not None)):
        raise Exception('either --libs_table or both --bam_files and --lib_ids are required')

    if (args['libs_table'] is not None) == (args['bam_files'] is not None or args['lib_ids'] is not None):
        raise Exception('--libs_table is mutually exclusive with --bam_files and --lib_ids')

    if args['bam_files'] is not None and (len(args['bam_files']) != len(args['lib_ids'])):
        raise Exception('--lib_ids must correspond one to one with --bam_files')

    config = {'ref_data_directory':args['ref_data_dir'],
              'package_data_directory':data_directory}
    execfile(default_config_filename, {}, config)

    if args['config'] is not None:
        execfile(args['config'], {}, config)

    config.update(args)

    pyp = pypeliner.app.Pypeline([destruct], config)


    # Pypeliner contexts

    locally = {'local':True}
    lowmem = {'mem':1}
    medmem = {'mem':8}
    himem = {'mem':16}

    workflow = pypeliner.workflow.Workflow()

    # Read in the bam file information

    if args['libs_table'] is not None:
        workflow.transform(
            name='readlibs',
            ctx=locally,
            func=destruct.read_libraries,
            ret=mgd.TempOutputObj('libinfo', 'bylibrary'),
            arg=(
                args['libs_table'],
            ),
        )
    else:
        workflow.transform(
            name='initlibs',
            ctx=locally,
            func=destruct.init_libraries,
            ret=mgd.TempOutputObj('libinfo', 'bylibrary'),
            args=(
                args['lib_ids'],
                args['bam_files'],
            ),
        )

    # Symlink the bam files locally

    workflow.transform(
        name='linklibs',
        axes=('bylibrary',),
        ctx=locally,
        func=destruct.link_libraries,
        args=(
            mgd.TempInputObj('libinfo', 'bylibrary').prop('bam'),
            mgd.TempOutputFile('bam', 'bylibrary'),
        ),
    )

    # Retrieve discordant reads and stats from bam files

    workflow.commandline(
        name='bamdisc',
        axes=('bylibrary',),
        ctx=medmem,
        args=(
            os.path.join(bin_directory, 'bamdiscordantfastq'),
            '-r',
            '-c', config['bam_max_soft_clipped'],
            '-f', config['bam_max_fragment_length'],
            '-b', mgd.TempInputFile('bam', 'bylibrary'),
            '-s', mgd.TempOutputFile('stats.file', 'bylibrary'),
            '-1', mgd.TempOutputFile('reads1', 'bylibrary'),
            '-2', mgd.TempOutputFile('reads2', 'bylibrary'),
            '-t', mgd.TempSpace('bamdisc.tempspace', 'bylibrary'),
        ),
    )

    workflow.commandline(
        name='bamsample',
        axes=('bylibrary',),
        ctx=medmem,
        args=(
            os.path.join(bin_directory, 'bamsamplefastq'),
            '-r',
            '-b', mgd.TempInputFile('bam', 'bylibrary'),
            '-n', config['num_read_samples'],
            '-1', mgd.TempOutputFile('sample1', 'bylibrary'),
            '-2', mgd.TempOutputFile('sample2', 'bylibrary'),
        ),
    )

    workflow.transform(
        name='readstats',
        axes=('bylibrary',),
        ctx=lowmem,
        func=destruct.read_stats,
        ret=mgd.TempOutputObj('stats', 'bylibrary'),
        args=(
            mgd.TempInputFile('stats.file', 'bylibrary'),
            config['fragment_length_num_stddevs'],
            mgd.TempOutputFile('flen.plots', 'bylibrary'),
            mgd.InputInstance('bylibrary'),
        ),
    )

    # Align a sample of reads and calculate alignment statistics

    workflow.transform(
        name='prepseed_sample',
        axes=('bylibrary',),
        ctx=medmem,
        func=destruct.prepare_seed_fastq,
        args=(
            mgd.TempInputFile('sample1', 'bylibrary'),
            mgd.TempInputFile('sample2', 'bylibrary'),
            36,
            mgd.TempOutputFile('sample.seed', 'bylibrary'),
        ),
    )

    workflow.commandline(
        name='bwtrealign_sample',
        axes=('bylibrary',),
        ctx=medmem,
        args=(
            'bowtie',
            config['genome_fasta'],
            mgd.TempInputFile('sample.seed', 'bylibrary'),
            '--chunkmbs', '512',
            '-k', '1000',
            '-m', '1000',
            '--strata',
            '--best',
            '-S',
            '|',
            os.path.join(bin_directory, 'aligntrue'),
            '-a', '-',
            '-1', mgd.TempInputFile('sample1', 'bylibrary'),
            '-2', mgd.TempInputFile('sample2', 'bylibrary'),
            '-r', config['genome_fasta'],
            '-g', config['gap_score'],
            '-x', config['mismatch_score'],
            '-m', config['match_score'],
            '--flmin', mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_min'),
            '--flmax', mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_max'),
            '-s', mgd.TempOutputFile('samples.align.true', 'bylibrary'),
        ),
    )

    workflow.transform(
        name='scorestats',
        axes=('bylibrary',),
        ctx=medmem,
        func=score_stats.create_score_stats,
        args=(
            mgd.TempInputFile('samples.align.true', 'bylibrary'),
            config['match_score'],
            mgd.TempOutputFile('score.stats', 'bylibrary'),
            mgd.TempOutputFile('score.stats.plots', 'bylibrary'),
            mgd.InputInstance('bylibrary'),
        ),
    )

    # Split discordant fastqs and align

    workflow.transform(
        name='splitfastq1',
        axes=('bylibrary',),
        ctx=lowmem,
        func=destruct.split_fastq,
        args=(
            mgd.TempInputFile('reads1', 'bylibrary'),
            int(config['reads_per_split']),
            mgd.TempOutputFile('reads1', 'bylibrary', 'byread'),
        ),
    )

    workflow.transform(
        name='splitfastq2',
        axes=('bylibrary',),
        ctx=lowmem,
        func=destruct.split_fastq,
        args=(
            mgd.TempInputFile('reads2', 'bylibrary'),
            int(config['reads_per_split']),
            mgd.TempOutputFile('reads2', 'bylibrary', 'byread', axes_origin=[]),
        ),
    )

    workflow.transform(
        name='prepseed',
        axes=('bylibrary', 'byread'),
        ctx=medmem,
        func=destruct.prepare_seed_fastq,
        args=(
            mgd.TempInputFile('reads1', 'bylibrary', 'byread'),
            mgd.TempInputFile('reads2', 'bylibrary', 'byread'),
            36,
            mgd.TempOutputFile('reads.seed', 'bylibrary', 'byread'),
        ),
    )

    workflow.commandline(
        name='bwtrealign',
        axes=('bylibrary', 'byread'),
        ctx=medmem,
        args=(
            'bowtie',
            config['genome_fasta'],
            mgd.TempInputFile('reads.seed', 'bylibrary', 'byread'),
            '--chunkmbs', '512',
            '-k', '1000',
            '-m', '1000',
            '--strata',
            '--best',
            '-S',
            '|',
            os.path.join(bin_directory, 'realign2'),
            '-l', mgd.TempInputObj('libinfo', 'bylibrary').prop('id'),
            '-a', '-',
            '-1', mgd.TempInputFile('reads1', 'bylibrary', 'byread'),
            '-2', mgd.TempInputFile('reads2', 'bylibrary', 'byread'),
            '-r', config['genome_fasta'],
            '-g', config['gap_score'],
            '-x', config['mismatch_score'],
            '-m', config['match_score'],
            '--flmin', mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_min'),
            '--flmax', mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_max'),
            '--tchimer', config['chimeric_threshold'],
            '--talign', config['alignment_threshold'],
            '--pchimer', config['chimeric_prior'],
            '--tvalid', config['readvalid_threshold'],
            '-z', mgd.TempInputFile('score.stats', 'bylibrary'),
            '--span', mgd.TempOutputFile('spanning.alignments', 'bylibrary', 'byread'),
            '--split', mgd.TempOutputFile('split.alignments', 'bylibrary', 'byread'),
        ),
    )

    workflow.transform(
        name='merge_spanning_1',
        axes=('bylibrary',),
        ctx=lowmem,
        func=destruct.merge_files_by_line,
        args=(
            mgd.TempInputFile('spanning.alignments', 'bylibrary', 'byread'),
            mgd.TempOutputFile('spanning.alignments_1', 'bylibrary'),
        ),
    )

    workflow.commandline(
        name='filterreads',
        axes=('bylibrary',),
        ctx=lowmem,
        args=(
            os.path.join(bin_directory, 'filterreads'),
            '-n', '2',
            '-a', mgd.TempInputFile('spanning.alignments_1', 'bylibrary'),
            '-r', config['satellite_regions'],
            '>', mgd.TempOutputFile('spanning.alignments', 'bylibrary'),
        ),
    )

    workflow.transform(
        name='merge_split_1',
        axes=('bylibrary',),
        ctx=lowmem,
        func=destruct.merge_files_by_line,
        args=(
            mgd.TempInputFile('split.alignments', 'bylibrary', 'byread'),
            mgd.TempOutputFile('split.alignments', 'bylibrary'),
        ),
    )

    workflow.transform(
        name='merge_spanning_2',
        ctx=lowmem,
        func=destruct.merge_files_by_line,
        args=(
            mgd.TempInputFile('spanning.alignments', 'bylibrary'),
            mgd.TempOutputFile('spanning.alignments'),
        ),
    )

    workflow.transform(
        name='merge_split_2',
        ctx=lowmem,
        func=destruct.merge_files_by_line,
        args=(
            mgd.TempInputFile('split.alignments', 'bylibrary'),
            mgd.TempOutputFile('split.alignments'),
        ),
    )

    # Cluster spanning reads

    workflow.transform(
        name='chromosome_args',
        ctx=locally,
        func=destruct.generate_chromosome_args,
        ret=mgd.TempOutputObj('chrom.args', 'bychromarg'),
        args=(
            config['chromosomes'],
        ),
    )

    workflow.transform(
        name='write_stats_table',
        ctx=lowmem,
        func=destruct.write_stats_table,
        args=(
            mgd.TempInputObj('libinfo', 'bylibrary'),
            mgd.TempInputObj('stats', 'bylibrary'),
            mgd.TempOutputFile('libstats.tsv'),
        ),
    )

    workflow.commandline(
        name='cluster',
        axes=('bychromarg',),
        ctx=medmem,
        args=(
            os.path.join(bin_directory, 'mclustermatepairs'),
            '-a', mgd.TempInputFile('spanning.alignments'),
            '-s', mgd.TempInputFile('libstats.tsv'),
            '-c', mgd.TempOutputFile('clusters', 'bychromarg'),
            mgd.TempInputObj('chrom.args', 'bychromarg'),
            '--clustmin', config['cluster_readcount_threshold'],
            '--fragmax', config['fragment_length_max'],
        ),
    )
    
    # Predict breakpoints from split reads

    workflow.transform(
        name='predict_breaks',
        axes=('bychromarg',),
        ctx=medmem,
        func=predict_breaks.predict_breaks,
        args=(
            mgd.TempInputFile('clusters', 'bychromarg'),
            mgd.TempInputFile('spanning.alignments'),
            mgd.TempInputFile('split.alignments'),
            mgd.TempOutputFile('breakpoints_2', 'bychromarg'),
        ),
    )

    workflow.transform(
        name='merge_clusters',
        ctx=lowmem,
        func=destruct.merge_clusters,
        args=(
            mgd.TempInputFile('clusters', 'bychromarg'),
            mgd.TempInputFile('breakpoints_2', 'bychromarg'),
            mgd.TempOutputFile('clusters'),
            mgd.TempOutputFile('breakpoints_2'),
            mgd.TempOutputFile('merge_clusters.debug'),
        ),
    )

    # Realign reads to breakpoints

    workflow.commandline(
        name='realigntobreaks',
        axes=('bylibrary', 'byread'),
        ctx=medmem,
        args=(
            os.path.join(bin_directory, 'realigntobreaks2'),
            '-r', config['genome_fasta'],
            '-b', mgd.TempInputFile('breakpoints_2'),
            '-c', mgd.TempInputFile('clusters'),
            '-g', config['gap_score'],
            '-x', config['mismatch_score'],
            '-m', config['match_score'],
            '--flmax', mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_max'),
            '--span', mgd.TempInputFile('spanning.alignments', 'bylibrary', 'byread'),
            '-1', mgd.TempInputFile('reads1', 'bylibrary', 'byread'),
            '-2', mgd.TempInputFile('reads2', 'bylibrary', 'byread'),
            '--realignments', mgd.TempOutputFile('realignments', 'bylibrary', 'byread'),
        ),
    )

    # Calculate likelihoods based on realignments

    workflow.transform(
        name='calculate_realignment_likelihoods',
        axes=('bylibrary', 'byread'),
        ctx=medmem,
        func=predict_breaks.calculate_realignment_likelihoods,
        args=(
            mgd.TempInputFile('breakpoints_2'),
            mgd.TempInputFile('realignments', 'bylibrary', 'byread'),
            mgd.TempInputFile('score.stats', 'bylibrary'),
            mgd.TempOutputFile('likelihoods_2', 'bylibrary', 'byread'),
            config['match_score'],
            mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_mean'),
            mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_stddev'),
        ),
    )

    workflow.transform(
        name='merge_likelihoods_1',
        axes=('bylibrary',),
        ctx=lowmem,
        func=destruct.merge_sorted_files_by_line,
        args=(
            mgd.TempInputFile('likelihoods_2', 'bylibrary', 'byread'),
            mgd.TempOutputFile('likelihoods_2', 'bylibrary'),
            mgd.TempSpace('merge_likelihoods_1_temp', 'bylibrary'),
            '1',
        ),
    )

    workflow.transform(
        name='merge_likelihoods_2',
        ctx=lowmem,
        func=destruct.merge_sorted_files_by_line,
        args=(
            mgd.TempInputFile('likelihoods_2', 'bylibrary'),
            mgd.TempOutputFile('likelihoods_2'),
            mgd.TempSpace('merge_likelihoods_2_temp'),
            '1',
        ),
    )

    # Set cover for multi mapping reads

    workflow.transform(
        name='calc_weights',
        ctx=medmem,
        func=predict_breaks.calculate_cluster_weights,
        args=(
            mgd.TempInputFile('breakpoints_2'),
            mgd.TempOutputFile('cluster_weights'),
        ),
    )

    workflow.commandline(
        name='setcover',
        ctx=medmem,
        args=(
            os.path.join(bin_directory, 'setcover'),
            '-c', mgd.TempInputFile('clusters'),
            '-w', mgd.TempInputFile('cluster_weights'),
            '-a', mgd.TempOutputFile('clusters_setcover'),
        ),
    )

    # Select cluster based on setcover

    workflow.transform(
        name='select_clusters',
        ctx=medmem,
        func=predict_breaks.select_clusters,
        args=(
            mgd.TempInputFile('clusters_setcover'),
            mgd.TempInputFile('breakpoints_2'),
            mgd.TempOutputFile('breakpoints_1'),
            mgd.TempInputFile('likelihoods_2'),
            mgd.TempOutputFile('likelihoods_1'),
        ),
    )

    # Select prediction based on max likelihood

    workflow.transform(
        name='select_predictions',
        ctx=himem,
        func=predict_breaks.select_predictions,
        args=(
            mgd.TempInputFile('breakpoints_1'),
            mgd.TempOutputFile('breakpoints'),
            mgd.TempInputFile('likelihoods_1'),
            mgd.TempOutputFile('likelihoods'),
            config['mate_score_threshold'],
            config['template_length_min_threshold'],
        ),
    )

    # Optionally tabulate supporting reads

    if args['breakpoint_read_table'] is not None:

        workflow.transform(
            name='tabreads',
            ctx=medmem,
            func=destruct.tabulate_reads,
            args=(
                mgd.TempInputFile('clusters_setcover'),
                mgd.TempInputObj('libinfo', 'bylibrary'),
                mgd.TempInputFile('reads1', 'bylibrary'),
                mgd.TempInputFile('reads2', 'bylibrary'),
                mgd.TempOutputFile('breakreads.table.unsorted'),
            ),
        )

        workflow.commandline(
            name='sortreads',
            ctx=medmem,
            args=(
                'sort', '-n',
                mgd.TempInputFile('breakreads.table.unsorted'),
                '>', mgd.OutputFile(args['breakpoint_read_table']),
            ),
        )


    # Tabulate results

    workflow.transform(
        name='tabulate',
        ctx=himem,
        func=destruct.tabulate_results,
        args=(
            mgd.TempInputFile('breakpoints'),
            mgd.TempInputFile('likelihoods'),
            mgd.TempInputObj('libinfo', 'bylibrary'),
            config['genome_fasta'],
            config['gtf_filename'],
            config['dgv_filename'],
            mgd.OutputFile(args['breakpoint_table']),
            mgd.OutputFile(args['breakpoint_library_table']),
        ),
    )

    workflow.transform(
        name='merge_plots',
        ctx=lowmem,
        func=destruct.merge_tars,
        args=(
            mgd.OutputFile(args['plots_tar']),
            mgd.TempInputFile('score.stats.plots', 'bylibrary'),
            mgd.TempInputFile('flen.plots', 'bylibrary'),
        ),
    )

    # Run the pipeline

    pyp.run(workflow)


else:


    def prepare_seed_fastq(reads_1_fastq, reads_2_fastq, seed_length, seed_fastq):
        with open(reads_1_fastq, 'r') as reads_1, open(reads_2_fastq, 'r') as reads_2, open(seed_fastq, 'w') as seed:
            fastq_lines = [[],[]]
            for fastq_1_line, fastq_2_line in zip(reads_1, reads_2):
                fastq_lines[0].append(fastq_1_line.rstrip())
                fastq_lines[1].append(fastq_2_line.rstrip())
                if len(fastq_lines[0]) == 4:
                    for read_end in (0,1):
                        if fastq_lines[read_end][1] > seed_length:
                            fastq_lines[read_end][1] = fastq_lines[read_end][1][0:seed_length]
                            fastq_lines[read_end][3] = fastq_lines[read_end][3][0:seed_length]
                        for line in fastq_lines[read_end]:
                            seed.write(line + '\n')
                    fastq_lines = [[],[]]


    class ConcordantReadStats(object):
        def __init__(self, stats, fragment_length_num_stddevs):
            self.stats = stats
            self.fragment_length_num_stddevs = fragment_length_num_stddevs
        @property
        def fragment_length_mean(self):
            return float(self.stats['fragment_mean'])
        @property
        def fragment_length_stddev(self):
            return float(self.stats['fragment_stddev'])
        @property
        def fragment_length_min(self):
            return int(self.fragment_length_mean - self.fragment_length_num_stddevs * self.fragment_length_stddev)
        @property
        def fragment_length_max(self):
            return int(self.fragment_length_mean + self.fragment_length_num_stddevs * self.fragment_length_stddev)


    def read_stats(stats_filename, fragment_length_num_stddevs, plots_tar_filename, library_id):
        stats = pd.read_csv(stats_filename, sep='\t')
        flen_stats = stats.loc[stats['type'] == 'fragment_length'].drop('type', axis=1)
        flen_stats = flen_stats.astype(float)
        fragment_count = flen_stats['value'].sum()
        fragment_mean = (flen_stats['key'] * flen_stats['value']).sum() / fragment_count
        fragment_variance = ((flen_stats['key'] - fragment_mean) * (flen_stats['key'] - fragment_mean) * flen_stats['value']).sum() / (fragment_count - 1)
        fragment_stddev = fragment_variance**0.5
        with tarfile.open(plots_tar_filename, 'w') as plots_tar:
            fig = plt.figure(figsize=(8,8))
            utils.plots.filled_density_weighted(plt.gca(), flen_stats['key'].values, flen_stats['value'].values, 'b', 0.5, 0, flen_stats['key'].max(), 4)
            plt.title('fragment lengths for library {0}'.format(library_id))
            utils.plots.savefig_tar(plots_tar, fig, 'fragment_length_{0}.pdf'.format(library_id))
            plt.clf()
        return ConcordantReadStats({'fragment_mean':fragment_mean, 'fragment_stddev':fragment_stddev}, fragment_length_num_stddevs)


    def write_stats_table(lib_infos, lib_stats, stats_table_filename):
        with open(stats_table_filename, 'w') as stats_table_file:
            for lib_name, lib_info in lib_infos.iteritems():
                stats_table_file.write(str(lib_infos[lib_name].id) + '\t')
                stats_table_file.write(str(lib_stats[lib_name].fragment_length_mean) + '\t')
                stats_table_file.write(str(lib_stats[lib_name].fragment_length_stddev) + '\n')


    def split_file_byline(in_filename, lines_per_file, out_filename_callback):
        with open(in_filename, 'r') as in_file:
            file_number = 0
            out_file = None
            out_file_lines = None
            try:
                for line in in_file:
                    if out_file is None or out_file_lines == lines_per_file:
                        if out_file is not None:
                            out_file.close()
                        out_file = open(out_filename_callback(file_number), 'w')
                        out_file_lines = 0
                        file_number += 1
                    out_file.write(line)
                    out_file_lines += 1
            finally:
                if out_file is not None:
                    out_file.close()


    def split_fastq(in_filename, num_reads_per_file, out_filename_callback):
        with open(in_filename, 'r') as in_file:
            file_number = 0
            out_file = None
            out_file_read_count = None
            try:
                for name, seq, comment, qual in itertools.izip_longest(*[in_file]*4):
                    if out_file is None or out_file_read_count == num_reads_per_file:
                        if out_file is not None:
                            out_file.close()
                        out_file = open(out_filename_callback(file_number), 'w')
                        out_file_read_count = 0
                        file_number += 1
                    out_file.write(name)
                    out_file.write(seq)
                    out_file.write(comment)
                    out_file.write(qual)
                    out_file_read_count += 1
            finally:
                if out_file is not None:
                    out_file.close()


    def merge_files_by_line(in_filenames, out_filename):
        with open(out_filename, 'w') as out_file:
            for id, in_filename in sorted(in_filenames.items()):
                with open(in_filename, 'r') as in_file:
                    for line in in_file:
                        out_file.write(line)


    def merge_sorted_files_by_line(in_filenames, out_filename, temp_space, sort_fields):
        pypeliner.commandline.execute(*['sort', '-T', temp_space, '-m', '-n', '-k', sort_fields] + in_filenames.values() + ['>', out_filename])


    def generate_chromosome_args(chromosomes):
        args = list()
        for chromosome_pair in itertools.combinations_with_replacement(chromosomes, 2):
            args.append('--inclchrompair ' + ','.join(chromosome_pair))
        args.append('--exclchrompairs ' + ','.join(chromosomes))
        return dict(enumerate(args))


    def read_clusters_breakpoints(clusters_filename, breakpoints_filename):
        with open(clusters_filename, 'r') as clusters_file, open(breakpoints_filename, 'r') as breakpoints_file:
            clusters_reader = csv.reader(clusters_file, delimiter='\t')
            breakpoints_reader = csv.reader(breakpoints_file, delimiter='\t')
            cluster_iter = itertools.groupby(clusters_reader, lambda row: row[0])
            breakend_iter = itertools.groupby(breakpoints_reader, lambda row: row[0])
            for (cluster_id_1, cluster_rows), (cluster_id_2, breakend_rows) in itertools.izip(cluster_iter, breakend_iter):
                if cluster_id_1 != cluster_id_2:
                    raise ValueError('Consistency issue between clusters and breakpoints for ' + clusters_filename + ' and ' + breakpoints_filename)
                yield cluster_id_1, cluster_rows, breakend_rows


    def merge_clusters(in_clusters_filenames, in_breakpoints_filenames,
                       out_clusters_filename, out_breakpoints_filename, debug_filename):
        new_cluster_id = 0
        with open(out_clusters_filename, 'w') as out_clusters_file, \
             open(out_breakpoints_filename, 'w') as out_breakpoints_file, \
             open(debug_filename, 'w') as debug_file:
            for idx, in_clusters_filename in in_clusters_filenames.iteritems():
                in_breakpoints_filename = in_breakpoints_filenames[idx]
                for cluster_id, cluster_rows, breakend_rows in read_clusters_breakpoints(in_clusters_filename, in_breakpoints_filename):
                    for row in cluster_rows:
                        row[0] = str(new_cluster_id)
                        out_clusters_file.write('\t'.join(row) + '\n')
                    for row in breakend_rows:
                        row[0] = str(new_cluster_id)
                        out_breakpoints_file.write('\t'.join(row) + '\n')
                    debug_file.write('{0}\t{1}\t{2}\n'.format(new_cluster_id, idx, cluster_id))
                    new_cluster_id += 1


    def remove_duplicates(input_filename, output_filename):
        with open(input_filename, 'r') as input_file, open(output_filename, 'w') as output_file:
            for cluster_id, cluster_rows in itertools.groupby(csv.reader(input_file, delimiter='\t'), lambda row: row[0]):
                fragment_mappings = set()
                for fragment_lib, fragment_rows in itertools.groupby(cluster_rows, lambda row: (row[2], row[11])):
                    fragment_rows = list(fragment_rows)
                    if len(fragment_rows) != 2:
                        raise Exception('require 2 lines per fragment for fragment ' + ','.join(fragment_lib))
                    fragment_mapping = set()
                    for read_end in (0, 1):
                        chromosome = fragment_rows[read_end][4]
                        position = fragment_rows[read_end][6] if fragment_rows[read_end][5] == '+' else fragment_rows[read_end][7]
                        fragment_mapping.add((fragment_lib[1], chromosome, position))
                    fragment_mapping = frozenset(fragment_mapping)
                    if fragment_mapping not in fragment_mappings:
                        for row in fragment_rows:
                            output_file.write('\t'.join(row) + '\n')
                        fragment_mappings.add(fragment_mapping)


    LibInfo = collections.namedtuple('LibInfo', ['id', 'name', 'bam'])


    def read_libraries(libraries_filename):
        libraries = dict()
        with open(libraries_filename, 'r') as libraries_file:
            for lib_idx, row in enumerate(csv.reader(libraries_file, delimiter='\t')):
                lib_name = row[0]
                lib_bam = row[1]
                libraries[lib_name] = LibInfo(lib_idx, lib_name, lib_bam)
        return libraries


    def init_libraries(lib_names, bam_filenames):
        libraries = dict()
        for lib_idx, (lib_name, lib_bam) in enumerate(sorted(zip(lib_names, bam_filenames))):
            libraries[lib_name] = LibInfo(lib_idx, lib_name, lib_bam)
        return libraries


    def link_libraries(target_bam_filename, link_bam_filename):
        try:
            os.remove(link_bam_filename)
        except OSError:
            pass
        os.symlink(os.path.abspath(target_bam_filename), link_bam_filename)


    def tabulate_reads(clusters_filename, lib_infos, reads1_filenames, reads2_filenames, reads_table_filename):
        fields = ['cluster_id', 'cluster_end', 'lib_id', 'read_id', 'read_end', 'align_id']
        clusters = pd.read_csv(clusters_filename, sep='\t', names=fields, usecols=['cluster_id', 'lib_id', 'read_id'])
        clusters = clusters.drop_duplicates().set_index(['lib_id', 'read_id']).sort_index()['cluster_id']
        with open(reads_table_filename, 'w') as reads_table_file:
            for lib_name in set(reads1_filenames.keys()).union(set(reads2_filenames.keys())):
                lib_id = lib_infos[lib_name].id
                for reads_filename in [reads1_filenames[lib_name], reads2_filenames[lib_name]]:
                    with open(reads_filename, 'r') as reads_file:
                        for name, seq, comment, qual in itertools.izip_longest(*[(a.rstrip() for a in reads_file)]*4):
                            assert name[0] == '@'
                            assert name[-1] == '1' or name[-1] == '2'
                            assert name[-2] == '/'
                            fragment_id = int(name[1:-2])
                            read_end = name[-1]
                            try:
                                cluster_ids = clusters.loc[(lib_id, fragment_id):(lib_id, fragment_id)]
                            except KeyError:
                                continue
                            for cluster_id in cluster_ids:
                                reads_table_file.write('\t'.join([str(cluster_id), str(lib_id), str(fragment_id), read_end, seq, qual, comment]) + '\n')


    class DGVDatabase(object):
        def __init__(self, dgv_filename):
            self.variations = list()
            chrvars = collections.defaultdict(list)
            with open(dgv_filename, 'r') as dgv_file:
                dgv_reader = csv.reader(dgv_file, delimiter='\t')
                dgv_header = next(dgv_reader)
                for row in dgv_reader:
                    id = row[0]
                    chr = row[1]
                    start = int(row[2])
                    end = int(row[3])
                    idx = len(self.variations)
                    chrvars[chr].append((idx, start, end))
                    self.variations.append((id, start, end))
            self.intervals = dict()
            for chr, vars in chrvars.iteritems():
                self.intervals[chr] = pygenes.IntervalTree(vars)
        def query(self, chromosome, start, end):
            if chromosome not in self.intervals:
                return
            idxs = self.intervals[chromosome].find_overlapping(start, end)
            for idx in [idxs[a] for a in range(0, len(idxs))]:
                startdiff = abs(start - self.variations[idx][1])
                enddiff = abs(end - self.variations[idx][2])
                if startdiff < 500 and enddiff < 500:
                    yield self.variations[idx][0]


    def merge_tars(output_filename, *input_filename_sets):
        with tarfile.open(output_filename, 'w') as output_tar:
            for input_filenames in input_filename_sets:
                for input_filename in input_filenames.itervalues():
                    with tarfile.open(input_filename, 'r') as in_tar:
                        for tarinfo in in_tar:
                            output_tar.addfile(tarinfo, in_tar.extractfile(tarinfo))


    converters = {'chromosome':str,
                  'chromosome_1':str,
                  'chromosome_2':str,
                  'inserted':str}


    def create_sequence(row, reference_sequences):
        breakend_sequences = ['', '']
        expected_strands = ('+', '-')
        inserted = ''
        if inserted != '.':
            inserted = row['inserted']
        for side in (0, 1):
            chromosome = row['chromosome_{0}'.format(side+1)]
            strand = row['strand_{0}'.format(side+1)]
            position = row['position_{0}'.format(side+1)]
            length = row['template_length_{0}'.format(side+1)]
            if strand == '+':
                start = position - length + 1
                end = position
            else:
                start = position
                end = position + length - 1
            breakend_sequences[side] = reference_sequences[chromosome][start-1:end]
            if strand != expected_strands[side]:
                breakend_sequences[side] = utils.misc.reverse_complement(breakend_sequences[side])
        return breakend_sequences[0] + '[' + inserted + ']' + breakend_sequences[1]


    def annotate_genes(row, gene_models):

        for side in (0, 1):

            chromosome = row['chromosome_{0}'.format(side+1)]
            position = row['position_{0}'.format(side+1)]

            nearest_gene_ids = gene_models.find_nearest_genes(chromosome, position)
            gene_id = 'NA'
            gene_name = 'NA'
            gene_location = 'NA'
            if len(nearest_gene_ids) > 0:
                gene_id = nearest_gene_ids[0]
                gene_name = gene_models.get_gene(gene_id).name
                gene_location = gene_models.calculate_gene_location(gene_id, position)

            row['gene_id_{0}'.format(side+1)] = gene_id
            row['gene_name_{0}'.format(side+1)] = gene_name
            row['gene_location_{0}'.format(side+1)] = gene_location

        return row


    def query_dgv(row, dgv):

        if row['chromosome_1'] != row['chromosome_2']:
            return 'NA'

        chromosome = row['chromosome_1']
        start, end = sorted((row['position_1'], row['position_2']))

        variants = list(dgv.query(chromosome, start, end))

        if len(variants) == 0:
            return 'NA'

        return ', '.join(variants)


    def tabulate_results(breakpoints_filename, likelihoods_filename, lib_infos,
                         genome_fasta, gtf_filename, dgv_filename,
                         breakpoint_table, breakpoint_library_table):

        lib_names = pd.DataFrame.from_records(lib_infos.values(), columns=LibInfo._fields, exclude=['bam'])
        lib_names = lib_names.rename(columns={'id':'library_id', 'name':'library'})

        breakpoints = pd.read_csv(breakpoints_filename, sep='\t',
                                  names=predict_breaks.breakpoint_fields,
                                  converters=converters)
        breakpoints = breakpoints.drop(['breakpoint_id'], axis=1)
        breakpoints = breakpoints.rename(columns={'count':'num_split'})
        breakpoints.loc[breakpoints['inserted'] == '.', 'inserted'] = ''

        breakpoints = utils.misc.normalize_breakpoints(breakpoints)

        likelihoods = pd.read_csv(likelihoods_filename, sep='\t',
                                  names=predict_breaks.likelihoods_fields,
                                  converters=converters)
        likelihoods = likelihoods.drop(['breakpoint_id'], axis=1)

        breakpoint_reads = (
            likelihoods.groupby(['cluster_id', 'library_id'])
            .size()
            .reset_index()
        )
        breakpoint_reads.columns = ['cluster_id', 'library_id', 'num_reads']

        breakpoint_unique_reads = (
            likelihoods.drop_duplicates(['cluster_id', 'library_id', 'template_length_1', 'template_length_2'])
            .groupby(['cluster_id', 'library_id'])
            .size()
            .reset_index()
        )
        breakpoint_unique_reads.columns = ['cluster_id', 'library_id', 'num_unique_reads']

        breakpoint_library = (
            breakpoint_reads.merge(breakpoint_unique_reads)
            .merge(lib_names)
            .drop(['library_id'], axis=1)
        )

        agg_f = {
            'log_likelihood':np.average,
            'log_cdf':np.average,
            'template_length_1':max,
            'template_length_2':max,
        }

        breakpoint_stats = (
            likelihoods.groupby('cluster_id')
            .agg(agg_f)
            .reset_index()
        )

        breakpoint_stats['template_length_min'] = breakpoint_stats[['template_length_1', 'template_length_2']].min(axis=1)

        breakpoint_counts = (
            likelihoods.groupby('cluster_id')
            .size()
            .reset_index()
        )
        breakpoint_counts.columns = ['cluster_id', 'num_reads']

        breakpoint_unique_counts = (
            likelihoods.drop_duplicates(['cluster_id', 'library_id', 'template_length_1', 'template_length_2'])
            .groupby('cluster_id')
            .size()
            .reset_index()
        )
        breakpoint_unique_counts.columns = ['cluster_id', 'num_unique_reads']

        breakpoints = breakpoints.merge(breakpoint_stats, on='cluster_id', how='left')

        breakpoints = breakpoints.merge(breakpoint_counts, on='cluster_id', how='left')

        # Calculate breakpoint type
        def breakpoint_type(row):
            if row['chromosome_1'] != row['chromosome_2']:
                return 'translocation'
            if row['strand_1'] == row['strand_2']:
                return 'inversion'
            positions = sorted([(row['position_{0}'.format(side)], row['strand_{0}'.format(side)]) for side in (1, 2)])
            if positions[0][1] == '+':
                return 'deletion'
            else:
                return 'duplication'

        breakpoints['type'] = breakpoints.apply(breakpoint_type, axis=1)

        # Calculate number inserted at the breakpoint
        def calculate_num_inserted(row):
            if row['inserted'] == '.':
                return 0
            else:
                return len(row['inserted'])

        breakpoints['num_inserted'] = breakpoints.apply(calculate_num_inserted, axis=1)

        # Annotate sequence
        reference_sequences = dict()
        for id, seq in utils.seq.read_sequences(open(genome_fasta, 'r')):
            reference_sequences[id] = seq

        breakpoints['sequence'] = breakpoints.apply(lambda row: create_sequence(row, reference_sequences), axis=1)

        # Annotate gene information
        gene_models = pygenes.GeneModels()
        gene_models.load_ensembl_gtf(gtf_filename)

        breakpoints = breakpoints.apply(lambda row: annotate_genes(row, gene_models), axis=1)

        # Annotate database of genomic variants
        dgv = DGVDatabase(dgv_filename)

        breakpoints['dgv_ids'] = breakpoints.apply(lambda row: query_dgv(row, dgv), axis=1)

        breakpoints = breakpoints.rename(columns={'cluster_id':'prediction_id'})

        breakpoints.to_csv(breakpoint_table, sep='\t', na_rep='NA', header=True, index=False)

        breakpoint_library = breakpoint_library.rename(columns={'cluster_id':'prediction_id'})

        breakpoint_library.to_csv(breakpoint_library_table, sep='\t', na_rep='NA', header=True, index=False)



