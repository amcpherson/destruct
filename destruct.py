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

import score_stats
import utils.plots
import utils.misc
import predict_breaks

__version__ = '0.1.0'

if __name__ == '__main__':

    import destruct

    argparser = argparse.ArgumentParser()
    pypeliner.easypypeliner.add_arguments(argparser)
    argparser.add_argument('--version', action='version', version=__version__)
    argparser.add_argument('libraries', help='Libraries list filename')
    argparser.add_argument('breakpoints', help='Breakpoints table filename')
    argparser.add_argument('breakreads', help='Breakpoint reads table filename')
    argparser.add_argument('plots_tar', help='Diagnostic plots tar filename')

    cfg = pypeliner.easypypeliner.Config(vars(argparser.parse_args()))
    pyp = pypeliner.easypypeliner.EasyPypeliner([destruct], cfg)

    pyp.sch.transform('readlibs', (), destruct.locally, destruct.read_libraries, pyp.sch.oobj('libinfo', ('bylibrary',)), cfg.libraries)

    pyp.sch.transform('linklibs', ('bylibrary',), destruct.locally, destruct.link_libraries, None, pyp.sch.iobj('libinfo', ('bylibrary',)).prop('bam'), pyp.sch.ofile('bam', ('bylibrary',)))

    destruct.multilib_predict_breakpoints(pyp.sch, cfg, pyp.sch.ifile('bam', ('bylibrary',)), pyp.sch.output(cfg.breakpoints), pyp.sch.output(cfg.breakreads), pyp.sch.output(cfg.plots_tar))

    pyp.run()

else:

    locally = {'local':True}
    lowmem = {'mem':1}
    medmem = {'mem':8}
    himem = {'mem':32}

    def multilib_predict_breakpoints(sch, cfg, bams, breakpoints, breakreads, plots_tar):

        '''
        inputs: 'bam', ('bylibrary',)
        outputs: 'breakpoints.table'
        '''

        retrieve_from_bam(sch, cfg, bams, ('bylibrary',))

        align_sample(sch, cfg, ('bylibrary',))
        split_align_merge(sch, cfg, ('bylibrary',))

        sch.transform('write_stats_table', (), lowmem,
            write_stats_table,
            None,
            sch.iobj('libinfo', ('bylibrary',)),
            sch.iobj('stats', ('bylibrary',)),
            sch.ofile('libstats.tsv'))

        sch.transform('merge_spanning', (), lowmem,
            merge_files_by_line,
            None,
            sch.ifile('spanning.alignments', ('bylibrary',)),
            sch.ofile('spanning.alignments'))

        sch.transform('merge_split', (), lowmem,
            merge_files_by_line,
            None,
            sch.ifile('split.alignments', ('bylibrary',)),
            sch.ofile('split.alignments'))


        # Cluster spanning reads

        sch.transform('chromosome_args', (), locally,
            generate_chromosome_args,
            sch.oobj('chrom.args', ('bychromarg',)),
            cfg.chromosomes.split(' '))

        sch.commandline('cluster', ('bychromarg',), medmem,
            cfg.mclustermatepairs_tool,
            '-a', sch.ifile('spanning.alignments'),
            '-s', sch.ifile('libstats.tsv'),
            '-c', sch.ofile('clusters', ('bychromarg',)),
            sch.iobj('chrom.args', ('bychromarg',)),
            '--clustmin', '1',
            '--fragmax', cfg.fragment_length_max)

        
        # Predict breakpoints from split reads

        sch.transform('predict_breaks', ('bychromarg',), medmem,
            predict_breaks.predict_breaks,
            None,
            sch.ifile('clusters', ('bychromarg',)),
            sch.ifile('spanning.alignments'),
            sch.ifile('split.alignments'),
            sch.ofile('breakpoints_2', ('bychromarg',)))

        sch.transform('merge_clusters', (), lowmem,
            merge_clusters,
            None,
            sch.ifile('clusters', ('bychromarg',)),
            sch.ifile('breakpoints_2', ('bychromarg',)),
            sch.ofile('clusters'),
            sch.ofile('breakpoints_2'),
            sch.ofile('merge_clusters.debug'))


        # Realign reads to breakpoints

        sch.commandline('realigntobreaks', ('bylibrary', 'byread'), lowmem,
            cfg.realigntobreaks2_tool,
            '-r', cfg.genome_fasta,
            '-b', sch.ifile('breakpoints_2'),
            '-c', sch.ifile('clusters'),
            '-g', cfg.gap_score,
            '-x', cfg.mismatch_score,
            '-m', cfg.match_score,
            '--flmax', sch.iobj('stats', ('bylibrary',)).prop('fragment_length_max'),
            '--span', sch.ifile('spanning.alignments', ('bylibrary', 'byread')),
            '--seqs', sch.ifile('reads', ('bylibrary', 'byread')),
            '--realignments', sch.ofile('realignments', ('bylibrary', 'byread')))


        # Calculate likelihoods based on realignments

        sch.transform('calculate_realignment_likelihoods', ('bylibrary', 'byread'), lowmem,
            predict_breaks.calculate_realignment_likelihoods,
            None,
            sch.ifile('breakpoints_2'),
            sch.ifile('realignments', ('bylibrary', 'byread')),
            sch.ifile('score.stats', ('bylibrary',)),
            sch.ofile('likelihoods_2', ('bylibrary', 'byread')),
            cfg.match_score,
            sch.iobj('stats', ('bylibrary',)).prop('fragment_length_mean'),
            sch.iobj('stats', ('bylibrary',)).prop('fragment_length_stddev'))

        sch.transform('merge_likelihoods1', ('bylibrary',), lowmem,
            merge_files_by_line,
            None,
            sch.ifile('likelihoods_2', ('bylibrary','byread')),
            sch.ofile('likelihoods_unsorted', ('bylibrary',)))

        sch.transform('merge_likelihoods2', (), lowmem,
            merge_files_by_line,
            None,
            sch.ifile('likelihoods_unsorted', ('bylibrary',)),
            sch.ofile('likelihoods_unsorted'))

        sch.commandline('sort_likelihoods', (), medmem,
            'sort', '-n', sch.ifile('likelihoods_unsorted'), '>', sch.ofile('likelihoods_2'))


        # Set cover for multi mapping reads

        sch.transform('calc_weights', (), lowmem,
            predict_breaks.calculate_cluster_weights,
            None,
            sch.ifile('breakpoints_2'),
            sch.ofile('cluster_weights'))

        sch.commandline('setcover', (), himem,
            cfg.setcover_tool,
            '-c', sch.ifile('clusters'),
            '-w', sch.ifile('cluster_weights'),
            '-a', sch.ofile('clusters_setcover'))


        # Select cluster based on setcover

        sch.transform('select_clusters', (), lowmem,
            predict_breaks.select_clusters,
            None,
            sch.ifile('clusters_setcover'),
            sch.ifile('breakpoints_2'),
            sch.ofile('breakpoints_1'),
            sch.ifile('likelihoods_2'),
            sch.ofile('likelihoods_1'))


        # Select prediction based on max likelihood

        sch.transform('select_predictions', (), lowmem,
            predict_breaks.select_predictions,
            None,
            sch.ifile('breakpoints_1'),
            sch.ofile('breakpoints'),
            sch.ifile('likelihoods_2'),
            sch.ofile('likelihoods'))


        # Predict rearrangement cycles

        sch.commandline('getclusterids', (), lowmem,
            'cut', '-f1', sch.ifile('clusters'), '|', 'uniq', '>', sch.ofile('clusters.ids'))

        sch.transform('splitclusterids', (), lowmem,
            split_file_byline,
            None,
            sch.ifile('clusters.ids'),
            int(cfg.clusters_per_split),
            sch.ofile('clusters.ids', ('bycluster',)))

        sch.commandline('cycles', ('bycluster',), himem,
            cfg.cycles_tool,
            '-b', sch.ifile('breakpoints'),
            '--idsfile', sch.ifile('clusters.ids', ('bycluster',)),
            '-s', cfg.cycles_scoremax,
            '-v', cfg.cycles_visitmax,
            '-y', cfg.cycles_lambda,
            '>', sch.ofile('cycles', ('bycluster',)))

        sch.transform('mergecycles', (), lowmem,
            merge_files_by_line,
            None,
            sch.ifile('cycles', ('bycluster',)),
            sch.ofile('cycles'))


        # Tabulate results

        sch.transform('tabreads', (), medmem,
            tabulate_reads,
            None,
            sch.ifile('clusters'),
            sch.iobj('libinfo', ('bylibrary',)),
            sch.ifile('reads1', ('bylibrary',)),
            sch.ifile('reads2', ('bylibrary',)),
            sch.ofile('breakreads.table.unsorted'))

        sch.commandline('sortreads', (), medmem,
            'sort', '-n', sch.ifile('breakreads.table.unsorted'), '>', breakreads)

        sch.transform('tabulate', (), himem,
            tabulate_results,
            None,
            sch.ifile('breakpoints'),
            sch.ifile('likelihoods'),
            sch.ifile('cycles'),
            sch.iobj('libinfo', ('bylibrary',)),
            cfg.genome_fasta,
            cfg.gtf_filename,
            cfg.dgv_filename,
            breakpoints)

        sch.transform('merge_plots', (), lowmem, merge_tars, None, plots_tar, sch.ifile('score.stats.plots', ('bylibrary',)), sch.ifile('flen.plots', ('bylibrary',)))


    def split_align_merge(sch, cfg, axes):

        '''
        inputs: 'reads1', 'reads2'
        outputs: 'spanning.alignments', 'split.alignments'
        '''

        sch.transform('splitfastq1', axes, lowmem, split_fastq, None, sch.ifile('reads1', axes), int(cfg.reads_per_split), sch.ofile('reads1', axes + ('byread',)))
        sch.transform('splitfastq2', axes, lowmem, split_fastq, None, sch.ifile('reads2', axes), int(cfg.reads_per_split), sch.ofile('reads2', axes + ('byread2',)))
        sch.changeaxis('fastq2axis', axes, 'reads2', 'byread2', 'byread')

        align_reads(sch, cfg, axes + ('byread',))

        sch.transform('mergespan', axes, lowmem, merge_files_by_line, None, sch.ifile('spanning.alignments', axes + ('byread',)), sch.ofile('spanning.alignments.unfiltered', axes))
        sch.commandline('filterreads', axes, lowmem, cfg.filterreads_tool, '-n', '2', '-a', sch.ifile('spanning.alignments.unfiltered', axes), '-r', cfg.satellite_regions, '>', sch.ofile('spanning.alignments', axes))
        sch.transform('mergesplt', axes, lowmem, merge_files_by_line, None, sch.ifile('split.alignments', axes + ('byread',)), sch.ofile('split.alignments', axes))


    def retrieve_from_bam(sch, cfg, bams, axes):

        '''
        inputs: 'bam'
        outputs: 'reads1', 'reads2', 'stats', 'sample1', 'sample2'
        '''

        sch.commandline('bamdisc', axes, medmem, cfg.bamdiscordantfastq_tool, '-r', '-c', cfg.bam_max_soft_clipped, '-f', cfg.bam_max_fragment_length, '-b', bams, '-s', sch.ofile('stats.file', axes), '-1', sch.ofile('reads1.unfiltered', axes), '-2', sch.ofile('reads2.unfiltered', axes), '-t', sch.tmpfile('bamdisc.tempspace', axes))
        sch.commandline('bamsample', axes, medmem, cfg.bamsamplefastq_tool, '-r', '-b', bams, '-n', cfg.num_read_samples, '-1', sch.ofile('sample1.unfiltered', axes), '-2', sch.ofile('sample2.unfiltered', axes))
        sch.transform('readstats', axes, lowmem, read_stats, sch.oobj('stats', axes), sch.ifile('stats.file', axes), float(cfg.fragment_length_num_stddevs), sch.ofile('flen.plots', axes), sch.inst('bylibrary'))
        sch.commandline('qtrimdisc', axes, lowmem, cfg.qualtrimfastq_tool, '-o', cfg.base_quality_offset, '-l', '36', '-q', '5', sch.ifile('reads1.unfiltered', axes), sch.ifile('reads2.unfiltered', axes), sch.ofile('reads1', axes), sch.ofile('reads2', axes))
        sch.commandline('qtrimsample', axes, lowmem, cfg.qualtrimfastq_tool, '-o', cfg.base_quality_offset, '-l', '36', '-q', '5', sch.ifile('sample1.unfiltered', axes), sch.ifile('sample2.unfiltered', axes), sch.ofile('sample1', axes), sch.ofile('sample2', axes))


    def align_sample(sch, cfg, axes):

        '''
        inputs: 'sample1', 'sample2'
        outputs: 'samples.align.true', 'samples.align.null'
        '''

        sch.commandline('bwtsample', axes, medmem, cfg.bowtie2_bin, '--no-mixed', '--no-discordant', '--very-sensitive', '-X', cfg.fragment_length_max, cfg.bowtie2_options, '-x', cfg.genome_fasta, '-1', sch.ifile('sample1', axes), '-2', sch.ifile('sample2', axes), '>', sch.ofile('sample.sam', axes))
        sch.commandline('aligntrue', axes, medmem, cfg.aligntrue_tool, '-g', cfg.gap_score, '-x', cfg.mismatch_score, '-m', cfg.match_score, '-r', cfg.genome_fasta, '-a', sch.ifile('sample.sam', axes), '>', sch.ofile('samples.align.true', axes))
        sch.transform('scorestats', axes, medmem, score_stats.create_score_stats, None, sch.ifile('samples.align.true', axes), int(cfg.match_score), sch.ofile('score.stats', axes), sch.ofile('score.stats.plots', axes), sch.inst('bylibrary'))


    def align_reads(sch, cfg, axes):

        '''
        inputs: 'reads1', 'reads2'
        outputs: 'spanning.alignments', 'split.alignments'
        '''

        sch.transform('prepseed', axes, medmem, prepare_seed_fastq, None, sch.ifile('reads1', axes), sch.ifile('reads2', axes), 36, sch.ofile('reads.seed', axes))
        sch.commandline('prepreads', axes, medmem, 'cat', sch.ifile('reads1', axes), sch.ifile('reads2', axes), '>', sch.ofile('reads', axes))
        sch.commandline('bwtrealign', axes, medmem, cfg.bowtie_bin, cfg.genome_fasta, sch.ifile('reads.seed', axes), '--chunkmbs', '512', '-k', '1000', '-m', '1000', '--strata', '--best', '-S', '|', cfg.realign2_tool, '-l', sch.iobj('libinfo', ('bylibrary',)).prop('id'), '-a', '-', '-s', sch.ifile('reads', axes), '-r', cfg.genome_fasta, '-g', cfg.gap_score, '-x', cfg.mismatch_score, '-m', cfg.match_score, '--flmin', sch.iobj('stats', ('bylibrary',)).prop('fragment_length_min'), '--flmax', sch.iobj('stats', ('bylibrary',)).prop('fragment_length_max'), '--tchimer', cfg.chimeric_threshold, '--talign', cfg.alignment_threshold, '--pchimer', cfg.chimeric_prior, '--tvalid', cfg.readvalid_threshold, '-z', sch.ifile('score.stats', ('bylibrary',)), '--span', sch.ofile('spanning.alignments', axes), '--split', sch.ofile('split.alignments', axes))


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


    def read_sequences(fasta):
        id = None
        sequences = []
        for line in fasta:
            line = line.rstrip()
            if len(line) == 0:
                continue
            if line[0] == '>':
                if id is not None:
                    yield (id, ''.join(sequences))
                id = line[1:]
                sequences = []
            else:
                sequences.append(line)
        if id is not None:
            yield (id, ''.join(sequences))


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
        breakpoints = breakpoints.drop(['prediction_id'], axis=1)
        breakpoints = breakpoints.rename(columns={'count':'num_split'})

        likelihoods = pd.read_csv(likelihoods_filename, sep='\t',
                                  names=predict_breaks.likelihoods_fields,
                                  converters=converters)
        likelihoods = likelihoods.drop(['prediction_id'], axis=1)

        agg_f = {'log_likelihood':sum,
                 'log_cdf':sum,
                 'template_length_1':max,
                 'template_length_2':max}

        breakpoint_stats = likelihoods.groupby('cluster_id').agg(agg_f)

        breakpoint_counts = likelihoods.groupby(['cluster_id', 'library_id'])\
                                       .size()\
                                       .unstack()\
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
        for id, seq in read_sequences(open(genome_fasta, 'r')):
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
                cycles_table.append((cluster_id, score_num_cycle_breaks, ', '.join(ids)))

        cycles_columns = ['cluster_id', 'cycle_score', 'cycle_num_breaks', 'cycle_ids']

        if len(cycles_table) == 0:
            cycles_table = pd.DataFrame(columns=cycles_columns)
        else:
            cycles_table = pd.DataFrame(cycles_table, columns=cycles_columns)

        breakpoints = breakpoints.merge(cycles_table, on='cluster_id', how='left')

        breakpoints.to_csv(results_filename, sep='\t', na_rep='NA')

