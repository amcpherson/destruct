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
import pypeliner.managed as mgd

import score_stats
import utils.plots
import utils.misc
import utils.io
import predict_breaks


destruct_directory = os.path.abspath(os.path.dirname(__file__))

data_directory = os.path.join(destruct_directory, 'data')
tools_directory = os.path.join(destruct_directory, 'tools')
default_config_filename = os.path.join(destruct_directory, 'defaultconfig.py')


__version__ = '0.1.0'

if __name__ == '__main__':

    import destruct

    argparser = argparse.ArgumentParser()

    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('--version', action='version', version=__version__)

    argparser.add_argument('ref_data_dir', help='Reference dataset directory')

    argparser.add_argument('breakpoint_table', help='Output table of breakpoint information in TSV format')

    argparser.add_argument('breakpoint_read_table', help='Output table of breakpoint read information in TSV format')

    argparser.add_argument('plots_tar', help='Output diagnostic plots tar filename')

    argparser.add_argument('--libs_table', help='Libraries list table filename')

    argparser.add_argument('--bam_files', nargs='+', help='Bam filenames')

    argparser.add_argument('--lib_ids', nargs='+', help='Ids for respective bam filenames')

    argparser.add_argument('--config', help='Configuration Filename')

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
    sch = pyp.sch


    # Pypeliner contexts

    locally = {'local':True}
    lowmem = {'mem':1}
    medmem = {'mem':8}
    himem = {'mem':16}


    # Read in the bam file information

    if args['libs_table'] is not None:
        sch.transform('readlibs', (), locally,
            destruct.read_libraries,
            mgd.TempOutputObj('libinfo', 'bylibrary'),
            args['libs_table'])
    else:
        sch.transform('initlibs', (), locally,
            destruct.init_libraries,
            mgd.TempOutputObj('libinfo', 'bylibrary'),
            args['lib_ids'],
            args['bam_files'])


    # Symlink the bam files locally

    sch.transform('linklibs', ('bylibrary',), locally,
        destruct.link_libraries,
        None,
        mgd.TempInputObj('libinfo', 'bylibrary').prop('bam'),
        mgd.TempOutputFile('bam', 'bylibrary'))


    # Retrieve discordant reads and stats from bam files

    sch.commandline('bamdisc', ('bylibrary',), medmem,
        os.path.join(tools_directory, 'bamdiscordantfastq'),
        '-r',
        '-c', config['bam_max_soft_clipped'],
        '-f', config['bam_max_fragment_length'],
        '-b', mgd.TempInputFile('bam', 'bylibrary'),
        '-s', mgd.TempOutputFile('stats.file', 'bylibrary'),
        '-1', mgd.TempOutputFile('reads1', 'bylibrary'),
        '-2', mgd.TempOutputFile('reads2', 'bylibrary'),
        '-t', mgd.TempFile('bamdisc.tempspace', 'bylibrary'))

    sch.commandline('bamsample', ('bylibrary',), medmem,
        os.path.join(tools_directory, 'bamsamplefastq'),
        '-r',
        '-b', mgd.TempInputFile('bam', 'bylibrary'),
        '-n', config['num_read_samples'],
        '-1', mgd.TempOutputFile('sample1', 'bylibrary'),
        '-2', mgd.TempOutputFile('sample2', 'bylibrary'))

    sch.transform('readstats', ('bylibrary',), lowmem,
        destruct.read_stats,
        mgd.TempOutputObj('stats', 'bylibrary'),
        mgd.TempInputFile('stats.file', 'bylibrary'),
        config['fragment_length_num_stddevs'],
        mgd.TempOutputFile('flen.plots', 'bylibrary'),
        mgd.InputInstance('bylibrary'))


    # Align a sample of reads and calculate alignment statistics

    sch.transform('prepseed_sample', ('bylibrary',), medmem, 
        destruct.prepare_seed_fastq,
        None,
        mgd.TempInputFile('sample1', 'bylibrary'),
        mgd.TempInputFile('sample2', 'bylibrary'),
        36,
        mgd.TempOutputFile('sample.seed', 'bylibrary'))

    sch.commandline('bwtrealign_sample', ('bylibrary',), medmem,
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
        os.path.join(tools_directory, 'aligntrue'),
        '-a', '-',
        '-1', mgd.TempInputFile('sample1', 'bylibrary'),
        '-2', mgd.TempInputFile('sample2', 'bylibrary'),
        '-r', config['genome_fasta'],
        '-g', config['gap_score'],
        '-x', config['mismatch_score'],
        '-m', config['match_score'],
        '--flmin', mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_min'),
        '--flmax', mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_max'),
        '-s', mgd.TempOutputFile('samples.align.true', 'bylibrary'))

    sch.transform('scorestats', ('bylibrary',), medmem,
        score_stats.create_score_stats,
        None,
        mgd.TempInputFile('samples.align.true', 'bylibrary'),
        config['match_score'],
        mgd.TempOutputFile('score.stats', 'bylibrary'),
        mgd.TempOutputFile('score.stats.plots', 'bylibrary'),
        mgd.InputInstance('bylibrary'))


    # Split discordant fastqs and align

    sch.transform('splitfastq1', ('bylibrary',), lowmem,
        destruct.split_fastq,
        None,
        mgd.TempInputFile('reads1', 'bylibrary'),
        int(config['reads_per_split']),
        mgd.TempOutputFile('reads1', 'bylibrary', 'byread'))

    sch.transform('splitfastq2', ('bylibrary',), lowmem,
        destruct.split_fastq,
        None,
        mgd.TempInputFile('reads2', 'bylibrary'),
        int(config['reads_per_split']),
        mgd.TempOutputFile('reads2', 'bylibrary', 'byread2'))

    sch.changeaxis('fastq2axis', ('bylibrary',), 'reads2', 'byread2', 'byread')

    sch.transform('prepseed', ('bylibrary', 'byread'), medmem, 
        destruct.prepare_seed_fastq,
        None,
        mgd.TempInputFile('reads1', 'bylibrary', 'byread'),
        mgd.TempInputFile('reads2', 'bylibrary', 'byread'),
        36,
        mgd.TempOutputFile('reads.seed', 'bylibrary', 'byread'))

    sch.commandline('bwtrealign', ('bylibrary', 'byread'), medmem,
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
        os.path.join(tools_directory, 'realign2'),
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
        '--split', mgd.TempOutputFile('split.alignments', 'bylibrary', 'byread'))

    sch.transform('mergespan', ('bylibrary',), lowmem,
        destruct.merge_files_by_line,
        None,
        mgd.TempInputFile('spanning.alignments', 'bylibrary', 'byread'),
        mgd.TempOutputFile('spanning.alignments_2', 'bylibrary'))

    sch.commandline('filterreads', ('bylibrary',), lowmem,
        os.path.join(tools_directory, 'filterreads'),
        '-n', '2',
        '-a', mgd.TempInputFile('spanning.alignments_2', 'bylibrary'),
        '-r', config['satellite_regions'],
        '>', mgd.TempOutputFile('spanning.alignments_1', 'bylibrary'))

    sch.transform('filterdups', ('bylibrary',), medmem,
        predict_breaks.remove_duplicates,
        None,
        mgd.TempInputFile('spanning.alignments_1', 'bylibrary'),
        mgd.TempOutputFile('spanning.alignments', 'bylibrary'))

    sch.transform('mergesplt', ('bylibrary',), lowmem,
        destruct.merge_files_by_line,
        None,
        mgd.TempInputFile('split.alignments', 'bylibrary', 'byread'),
        mgd.TempOutputFile('split.alignments', 'bylibrary'))

    sch.transform('merge_spanning', (), lowmem,
        destruct.merge_files_by_line,
        None,
        mgd.TempInputFile('spanning.alignments', 'bylibrary'),
        mgd.TempOutputFile('spanning.alignments'))

    sch.transform('merge_split', (), lowmem,
        destruct.merge_files_by_line,
        None,
        mgd.TempInputFile('split.alignments', 'bylibrary'),
        mgd.TempOutputFile('split.alignments'))


    # Cluster spanning reads

    sch.transform('chromosome_args', (), locally,
        destruct.generate_chromosome_args,
        mgd.TempOutputObj('chrom.args', 'bychromarg'),
        config['chromosomes'])

    sch.transform('write_stats_table', (), lowmem,
        destruct.write_stats_table,
        None,
        mgd.TempInputObj('libinfo', 'bylibrary'),
        mgd.TempInputObj('stats', 'bylibrary'),
        mgd.TempOutputFile('libstats.tsv'))

    sch.commandline('cluster', ('bychromarg',), medmem,
        os.path.join(tools_directory, 'mclustermatepairs'),
        '-a', mgd.TempInputFile('spanning.alignments'),
        '-s', mgd.TempInputFile('libstats.tsv'),
        '-c', mgd.TempOutputFile('clusters', 'bychromarg'),
        mgd.TempInputObj('chrom.args', 'bychromarg'),
        '--clustmin', config['cluster_readcount_threshold'],
        '--fragmax', config['fragment_length_max'])

    
    # Predict breakpoints from split reads

    sch.transform('predict_breaks', ('bychromarg',), medmem,
        predict_breaks.predict_breaks,
        None,
        mgd.TempInputFile('clusters', 'bychromarg'),
        mgd.TempInputFile('spanning.alignments'),
        mgd.TempInputFile('split.alignments'),
        mgd.TempOutputFile('breakpoints_2', 'bychromarg'))

    sch.transform('merge_clusters', (), lowmem,
        destruct.merge_clusters,
        None,
        mgd.TempInputFile('clusters', 'bychromarg'),
        mgd.TempInputFile('breakpoints_2', 'bychromarg'),
        mgd.TempOutputFile('clusters'),
        mgd.TempOutputFile('breakpoints_2'),
        mgd.TempOutputFile('merge_clusters.debug'))


    # Realign reads to breakpoints

    sch.commandline('realigntobreaks', ('bylibrary', 'byread'), medmem,
        os.path.join(tools_directory, 'realigntobreaks2'),
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
        '--realignments', mgd.TempOutputFile('realignments', 'bylibrary', 'byread'))


    # Calculate likelihoods based on realignments

    sch.transform('calculate_realignment_likelihoods', ('bylibrary', 'byread'), medmem,
        predict_breaks.calculate_realignment_likelihoods,
        None,
        mgd.TempInputFile('breakpoints_2'),
        mgd.TempInputFile('realignments', 'bylibrary', 'byread'),
        mgd.TempInputFile('score.stats', 'bylibrary'),
        mgd.TempOutputFile('likelihoods_2', 'bylibrary', 'byread'),
        config['match_score'],
        mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_mean'),
        mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_stddev'))

    sch.transform('merge_likelihoods1', ('bylibrary',), lowmem,
        destruct.merge_files_by_line,
        None,
        mgd.TempInputFile('likelihoods_2', 'bylibrary', 'byread'),
        mgd.TempOutputFile('likelihoods_2', 'bylibrary'))

    sch.transform('merge_likelihoods2', (), lowmem,
        destruct.merge_files_by_line,
        None,
        mgd.TempInputFile('likelihoods_2', 'bylibrary'),
        mgd.TempOutputFile('likelihoods_2'))


    # Set cover for multi mapping reads

    sch.transform('calc_weights', (), medmem,
        predict_breaks.calculate_cluster_weights,
        None,
        mgd.TempInputFile('breakpoints_2'),
        mgd.TempOutputFile('cluster_weights'))

    sch.commandline('setcover', (), medmem,
        os.path.join(tools_directory, 'setcover'),
        '-c', mgd.TempInputFile('clusters'),
        '-w', mgd.TempInputFile('cluster_weights'),
        '-a', mgd.TempOutputFile('clusters_setcover'))


    # Select cluster based on setcover

    sch.transform('select_clusters', (), medmem,
        predict_breaks.select_clusters,
        None,
        mgd.TempInputFile('clusters_setcover'),
        mgd.TempInputFile('breakpoints_2'),
        mgd.TempOutputFile('breakpoints_1'),
        mgd.TempInputFile('likelihoods_2'),
        mgd.TempOutputFile('likelihoods_1'))


    # Select prediction based on max likelihood

    sch.transform('select_predictions', (), himem,
        predict_breaks.select_predictions,
        None,
        mgd.TempInputFile('breakpoints_1'),
        mgd.TempOutputFile('breakpoints'),
        mgd.TempInputFile('likelihoods_1'),
        mgd.TempOutputFile('likelihoods'),
        config['mate_score_threshold'],
        config['template_length_min_threshold'])
    

    # Predict rearrangement cycles

    sch.commandline('getclusterids', (), lowmem,
        'cut', '-f1', mgd.TempInputFile('clusters'), '|', 'uniq', '>', mgd.TempOutputFile('clusters.ids'))

    sch.transform('splitclusterids', (), lowmem,
        destruct.split_file_byline,
        None,
        mgd.TempInputFile('clusters.ids'),
        config['clusters_per_split'],
        mgd.TempOutputFile('clusters.ids', 'bycluster'))

    sch.commandline('cycles', ('bycluster',), medmem,
        os.path.join(tools_directory, 'cycles'),
        '-b', mgd.TempInputFile('breakpoints'),
        '--idsfile', mgd.TempInputFile('clusters.ids', 'bycluster'),
        '-s', config['cycles_scoremax'],
        '-v', config['cycles_visitmax'],
        '-y', config['cycles_lambda'],
        '>', mgd.TempOutputFile('cycles', 'bycluster'))

    sch.transform('mergecycles', (), lowmem,
        destruct.merge_files_by_line,
        None,
        mgd.TempInputFile('cycles', 'bycluster'),
        mgd.TempOutputFile('cycles'))


    # Tabulate results

    sch.transform('tabreads', (), medmem,
        destruct.tabulate_reads,
        None,
        mgd.TempInputFile('clusters_setcover'),
        mgd.TempInputObj('libinfo', 'bylibrary'),
        mgd.TempInputFile('reads1', 'bylibrary'),
        mgd.TempInputFile('reads2', 'bylibrary'),
        mgd.TempOutputFile('breakreads.table.unsorted'))

    sch.commandline('sortreads', (), medmem,
        'sort', '-n', mgd.TempInputFile('breakreads.table.unsorted'), '>', mgd.OutputFile(args['breakpoint_read_table']))

    sch.transform('tabulate', (), himem,
        destruct.tabulate_results,
        None,
        mgd.TempInputFile('breakpoints'),
        mgd.TempInputFile('likelihoods'),
        mgd.TempInputFile('cycles'),
        mgd.TempInputObj('libinfo', 'bylibrary'),
        config['genome_fasta'],
        config['gtf_filename'],
        config['dgv_filename'],
        mgd.OutputFile(args['breakpoint_table']))

    sch.transform('merge_plots', (), lowmem,
        destruct.merge_tars,
        None,
        mgd.OutputFile(args['plots_tar']),
        mgd.TempInputFile('score.stats.plots', 'bylibrary'),
        mgd.TempInputFile('flen.plots', 'bylibrary'))


    # Run the pipeline

    pyp.run()


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
                  'chromosome_2':str}


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


    def tabulate_results(breakpoints_filename, likelihoods_filename,
                         cycles_filename, lib_infos,
                         genome_fasta, gtf_filename, dgv_filename,
                         results_filename):

        lib_names = dict([(info.id, info.name) for info in lib_infos.values()])

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

        agg_f = {'log_likelihood':np.average,
                 'log_cdf':np.average,
                 'template_length_1':max,
                 'template_length_2':max}

        breakpoint_stats = likelihoods.groupby('cluster_id').agg(agg_f)

        breakpoint_stats['template_length_min'] = breakpoint_stats[['template_length_1', 'template_length_2']].min(axis=1)

        breakpoint_counts = likelihoods.groupby(['cluster_id', 'library_id'])\
                                       .size()\
                                       .unstack()\
                                       .fillna(0)\
                                       .astype(int)\
                                       .rename(columns=lib_names)
        breakpoint_counts.columns = [a+'_count' for a in breakpoint_counts.columns]

        breakpoints = breakpoints.merge(breakpoint_stats, left_on='cluster_id', right_index=True, how='left')

        breakpoints = breakpoints.merge(breakpoint_counts, left_on='cluster_id', right_index=True, how='left')

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
        for id, seq in utils.io.read_sequences(open(genome_fasta, 'r')):
            reference_sequences[id] = seq

        breakpoints['sequence'] = breakpoints.apply(lambda row: create_sequence(row, reference_sequences), axis=1)

        # Annotate gene information
        gene_models = pygenes.GeneModels()
        gene_models.load_ensembl_gtf(gtf_filename)

        breakpoints = breakpoints.apply(lambda row: annotate_genes(row, gene_models), axis=1)

        # Annotate database of genomic variants
        dgv = DGVDatabase(dgv_filename)

        breakpoints['dgv_ids'] = breakpoints.apply(lambda row: query_dgv(row, dgv), axis=1)

        # Annotate rearrangement cycles
        cycles_table = list()
        with open(cycles_filename, 'r') as cycles_file:
            for row in csv.reader(cycles_file, delimiter='\t'):
                score = row[0]
                ids = row[1::4]
                cluster_id = ids[0]
                num_cycle_breaks = len(ids)
                cycles_table.append((cluster_id, score, num_cycle_breaks, ', '.join(ids)))

        cycles_columns = ['cluster_id', 'cycle_score', 'cycle_num_breaks', 'cycle_ids']

        if len(cycles_table) == 0:
            cycles_table = pd.DataFrame(columns=cycles_columns)
        else:
            cycles_table = pd.DataFrame(cycles_table, columns=cycles_columns)

        breakpoints = breakpoints.merge(cycles_table, on='cluster_id', how='left')

        breakpoints = breakpoints.rename(columns={'cluster_id':'prediction_id'})

        breakpoints.to_csv(results_filename, sep='\t', na_rep='NA', header=True, index=False)


