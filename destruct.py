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
import pandas as pd
import matplotlib.pyplot as plt

import pygenes
import pypeliner

import score_stats
import utils.plots

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

        sch.transform('write_stats_table', (), lowmem, write_stats_table, None, sch.iobj('libinfo', ('bylibrary',)), sch.iobj('stats', ('bylibrary',)), sch.ofile('libstats.tsv'))

        sch.transform('merge_alignments', (), lowmem, merge_files_by_line, None, sch.ifile('spanning.alignments', ('bylibrary',)), sch.ofile('spanning.alignments'))

        sch.transform('chromosome_args', (), locally, generate_chromosome_args, sch.oobj('chrom.args', ('bychromarg',)), cfg.chromosomes.split(' '))

        sch.commandline('cluster', ('bychromarg',), medmem, cfg.mclustermatepairs_tool, '-a', sch.ifile('spanning.alignments'), '-s', sch.ifile('libstats.tsv'), '-c', sch.ofile('clusters.raw', ('bychromarg',)), sch.iobj('chrom.args', ('bychromarg',)), '--clustmin', '1', '--fragmax', cfg.fragment_length_max)

        sch.transform('merge_clusters', (), lowmem, merge_clusters, None, sch.ifile('clusters.raw', ('bychromarg',)), sch.ofile('clusters.raw'), sch.ofile('merge_clusters.debug'))

        sch.transform('breaks', (), himem, run_mpredictbreaks, None, cfg.mpredictbreaks_tool, cfg.genome_fasta, sch.ifile('clusters.raw'), sch.ifile('split.alignments', ('bylibrary',)), sch.ofile('breakpoints'))

        sch.commandline('rankclust', (), himem, cfg.rankclusters_tool, '-c', sch.ifile('clusters.nodup'), '>', sch.ofile('clusters.prob'))

        sch.transform('calc_weights', (), lowmem, calculate_cluster_weights, None, sch.ifile('breakpoints'), sch.ofile('clusters.weights'))
        sch.commandline('setcover', (), himem, cfg.setcover_tool, '-c', sch.ifile('clusters.nodup'), '-w', sch.ifile('clusters.weights'), '-a', sch.ofile('clusters.setcover'))

        sch.transform('filter', (), medmem, filter_clusters, None, sch.ofile('clusters.filtered'), sch.ifile('clusters.setcover'), sch.ifile('clusters.prob'), float(cfg.cluster_align_threshold), float(cfg.cluster_chimeric_threshold), float(cfg.cluster_valid_threshold), int(cfg.cluster_coverage_threshold), int(cfg.cluster_readcount_threshold))

        sch.commandline('getclusterids', (), lowmem, 'cut', '-f1', sch.ifile('clusters.filtered'), '|', 'uniq', '>', sch.ofile('clusters.filtered.ids'))
        sch.transform('splitclusterids', (), lowmem, split_file_byline, None, sch.ifile('clusters.filtered.ids'), int(cfg.clusters_per_split), sch.ofile('clusters.filtered.ids', ('bycluster',)))
        sch.commandline('cycles', ('bycluster',), himem, cfg.cycles_tool, '-c', sch.ifile('clusters.filtered'), '-p', sch.ifile('clusters.prob'), '--idsfile', sch.ifile('clusters.filtered.ids', ('bycluster',)), '-s', cfg.cycles_scoremax, '-v', cfg.cycles_visitmax, '-y', cfg.cycles_lambda, '>', sch.ofile('cycles', ('bycluster',)))
        sch.transform('mergecycles', (), lowmem, merge_files_by_line, None, sch.ifile('cycles', ('bycluster',)), sch.ofile('cycles'))

        sch.transform('tabreads', (), medmem, tabulate_reads, None, sch.ifile('clusters.filtered'), sch.ifile('reads1', ('bylibrary',)), sch.ifile('reads2', ('bylibrary',)), sch.ofile('breakreads.table.unsorted'))
        sch.commandline('sortreads', (), medmem, 'sort', '-n', sch.ifile('breakreads.table.unsorted'), '>', breakreads)

        sch.transform('tabulate', (), himem, multilib_tabulate, None, breakpoints, sch.ifile('clusters.filtered'), sch.ifile('clusters.prob'), sch.ifile('clusters.nodup'), sch.ifile('breakpoints'), cfg.genome_fasta, cfg.gtf_filename, cfg.dgv_filename, sch.ifile('cycles'), sch.iobj('stats', ('bylibrary',)))
        
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
            args.append('--chrompair ' + ','.join(chromosome_pair))
        args.append('--exclchrompairs ' + ','.join(chromosomes))
        return dict(enumerate(args))


    def merge_clusters(in_filenames, out_filename, debug_filename):
        new_cluster_id = 0
        with open(out_filename, 'w') as out_file, open(debug_filename, 'w') as debug_file:
            for idx, in_filename in in_filenames.iteritems():
                with open(in_filename, 'r') as in_file:
                    reader = csv.reader(in_file, delimiter='\t')
                    for cluster_id, rows in itertools.groupby(reader, lambda row: row[0]):
                        for row in rows:
                            row[0] = str(new_cluster_id)
                            out_file.write('\t'.join(row) + '\n')
                        debug_file.write('{0}\t{1}\t{2}\n'.format(new_cluster_id, idx, cluster_id))
                        new_cluster_id += 1


    def calculate_cluster_weights(breakpoints_filename, weights_filename):
        epsilon = 0.0001
        itx_distance = 1000000000
        with open(breakpoints_filename, 'r') as breakpoints_file, open(weights_filename, 'w') as weights_file:
            for row in csv.reader(breakpoints_file, delimiter='\t'):
                cluster_id = row[0]
                chromosome1 = row[1]
                chromosome2 = row[4]
                position1 = int(row[3])
                position2 = int(row[6])
                if chromosome1 != chromosome2:
                    distance = itx_distance
                else:
                    distance = abs(position1 - position2)
                weight = epsilon * math.log(distance)
                weights_file.write('\t'.join((cluster_id, str(weight))) + '\n')


    def merge_samples(all_samples, merged_samples):
        samples = []
        for sample_filename in all_samples.values():
            with open(sample_filename, 'r') as sample_file:
                for idx, line in enumerate(sample_file):
                    if idx >= len(samples):
                        samples.extend([0] * (idx - len(samples) + 1))
                    samples[idx] += int(line)
        with open(merged_samples, 'w') as merged_file:
            for sample in samples:
                merged_file.write('%d\n' % (sample,))


    def segregate_mitochondrial(mitochondrial_chromosome, input_filename, output_filename):
        with open(input_filename, 'r') as input_file, open(output_filename, 'w') as output_file:
            for cluster_id, rows in itertools.groupby(csv.reader(input_file, delimiter='\t'), lambda row: row[0]):
                rows = list(rows)
                chromosomes = [None,None]
                for row in rows:
                    chromosomes[int(row[1])] = row[4]
                if chromosomes[0] != chromosomes[1] and (chromosomes[0] == mitochondrial_chromosome or chromosomes[1] == mitochondrial_chromosome):
                    continue
                for row in rows:
                    output_file.write('\t'.join(row) + '\n')


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


    def run_mpredictbreaks(mpredictbreaks_bin, genome_fasta, clusters_filename, split_alignments, breakpoints):
        mpredictbreaks_arguments = [mpredictbreaks_bin, '-r', genome_fasta, '-c', clusters_filename, '-b', breakpoints, '-l', '-']
        mpredictbreaks_process = subprocess.Popen(mpredictbreaks_arguments, stdin=subprocess.PIPE, stdout=sys.stdout, stderr=sys.stderr)
        for lib_name, lib_split_alignments in split_alignments.iteritems():
            mpredictbreaks_process.stdin.write('\t'.join([lib_name, lib_split_alignments]) + '\n')
        mpredictbreaks_process.stdin.close()
        if mpredictbreaks_process.wait() != 0:
            raise Exception('mpredictbreaks process %s produced exit code %s' % (' '.join(mpredictbreaks_arguments), mpredictbreaks_process.returncode))


    def read_cluster_fragments(input_file):
        for row in csv.reader(input_file, delimiter='\t'):
            cluster_id = int(row[0])
            fragment_id = int(row[2])
            lib_id = row[11]
            yield cluster_id, fragment_id, lib_id


    def tabulate_reads(clusters_filename, reads1_filenames, reads2_filenames, reads_table_filename):
        cluster_fragments = defaultdict(set)
        with open(clusters_filename, 'r') as clusters_file:
            for cluster_id, fragment_id, lib_id in read_cluster_fragments(clusters_file):
                cluster_fragments[(fragment_id, lib_id)].add(cluster_id)
        with open(reads_table_filename, 'w') as reads_table_file:
            for lib_id in set(reads1_filenames.keys()).union(set(reads2_filenames.keys())):
                for reads_filename in [reads1_filenames[lib_id], reads2_filenames[lib_id]]:
                    with open(reads_filename, 'r') as reads_file:
                        for name, seq, comment, qual in itertools.izip_longest(*[(a.rstrip() for a in reads_file)]*4):
                            assert name[0] == '@'
                            assert name[-1] == '1' or name[-1] == '2'
                            assert name[-2] == '/'
                            fragment_id = int(name[1:-2])
                            read_end = name[-1]
                            for cluster_id in cluster_fragments.get((fragment_id, lib_id), set([])):
                                reads_table_file.write('\t'.join([str(cluster_id), lib_id, str(fragment_id), read_end, seq, qual, comment]) + '\n')


    GenomicRegion = collections.namedtuple('GenomicRegion', ['chromosome', 'strand', 'start', 'end'])

    def read_cluster_regions(input_file):
        for cluster_id, rows in itertools.groupby(csv.reader(input_file, delimiter='\t'), lambda row: row[0]):
            rows = list(rows)
            chromosomes = [None, None]
            strands = [None, None]
            starts = [[], []]
            ends = [[], []]
            lib_counts = Counter()
            for row in rows:
                cluster_end = int(row[1])
                chromosomes[cluster_end] = row[4]
                strands[cluster_end] = row[5]
                starts[cluster_end].append(int(row[6]))
                ends[cluster_end].append(int(row[7]))
                if cluster_end == 0:
                    lib_counts[row[11]] += 1
            regions = list()
            for cluster_end in (0, 1):
                regions.append(GenomicRegion(chromosomes[cluster_end], strands[cluster_end], min(starts[cluster_end]), max(ends[cluster_end])))
            yield cluster_id, regions, lib_counts, rows


    def read_breakpoints(input_file):
        for cluster_id, rows in itertools.groupby(csv.reader(input_file, delimiter='\t'), lambda row: row[0]):
            rows = list(rows)
            homology = 0
            prev_pos = None
            prev_dist = None
            for row in sorted(rows, key=lambda row: int(row[3])):
                current_pos = int(row[3])
                if row[2] == row[5]:
                    current_dist = int(row[3]) + int(row[6])
                else:
                    current_dist = int(row[3]) - int(row[6])
                if prev_pos is not None:
                    if current_pos == prev_pos + 1 and current_dist == prev_dist:
                        homology += 1
                    else:
                        homology = 0
                prev_pos = current_pos
                prev_dist = current_dist
            yield cluster_id, [int(rows[0][3]), int(rows[0][6])], int(rows[0][7]), int(rows[0][8]), homology


    def clusters_table_header(lib_names):
        header = 'cluster_id\t'
        for cluster_end in ('1', '2'):
            header += 'chromosome_' + cluster_end + '\t'
            header += 'strand_' + cluster_end + '\t'
            header += 'start_' + cluster_end + '\t'
            header += 'end_' + cluster_end + '\t'
            header += 'gene_id' + cluster_end + '\t'
            header += 'gene_name' + cluster_end + '\t'
            header += 'gene_location' + cluster_end + '\t'
            header += 'break_' + cluster_end + '\t'
        for lib_name in lib_names:
            header += lib_name + '_count' + '\t'
        header += 'exact_break' + '\t'
        header += 'align_prob' + '\t'
        header += 'chimeric_prob' + '\t'
        header += 'valid_prob' + '\t'
        header += 'type' + '\t'
        header += 'setcover_ratio' + '\t'
        header += 'num_split' + '\t'
        header += 'num_inserted' + '\t'
        header += 'homology' + '\t'
        header += 'cycle_score' + '\t'
        header += 'cycle_num_breaks' + '\t'
        header += 'cycle_ids' + '\t'
        header += 'dgv_ids' + '\t'
        header += 'sequence_approx' + '\t'
        header += 'sequence_exact'
        header += '\n'
        return header


    def clusters_table_row(lib_names, cluster_id, cluster_info, probs, setcover_ratio, breakpoint_info, dgv_info, cycle_info, sequences, gene_models):
        exact_break = 1
        break_positions = breakpoint_info[3]
        if break_positions is None:
            exact_break = 0
            break_positions = []
            for cluster_end in (0, 1):
                break_positions.append((cluster_info[0][cluster_end].start, cluster_info[0][cluster_end].end)[cluster_info[0][cluster_end].strand == '+'])
        row = cluster_id + '\t'
        for cluster_end in (0, 1):
            row += cluster_info[0][cluster_end].chromosome + '\t'
            row += cluster_info[0][cluster_end].strand + '\t'
            row += str(cluster_info[0][cluster_end].start) + '\t'
            row += str(cluster_info[0][cluster_end].end) + '\t'
            mid = (cluster_info[0][cluster_end].start + cluster_info[0][cluster_end].end) / 2
            nearest_gene_ids = gene_models.find_nearest_genes(cluster_info[0][cluster_end].chromosome, mid)
            gene_id = 'NA'
            gene_name = 'NA'
            gene_location = 'NA'
            if len(nearest_gene_ids) > 0:
                gene_id = nearest_gene_ids[0]
                gene_name = gene_models.get_gene(gene_id).name
                gene_location = gene_models.calculate_gene_location(gene_id, mid)
            row += gene_id + '\t'
            row += gene_name + '\t'
            row += gene_location + '\t'
            row += str(break_positions[cluster_end]) + '\t'
        for lib_name in lib_names:
            row += str(cluster_info[1].get(lib_name, 0)) + '\t'
        row += str(exact_break) + '\t'
        row += str(probs[cluster_id][0]) + '\t'
        row += str(probs[cluster_id][1]) + '\t'
        row += str(probs[cluster_id][2]) + '\t'
        row += rearrangement_type(cluster_info[0]) + '\t'
        row += str(setcover_ratio) + '\t'
        row += str(breakpoint_info[0]) + '\t'
        row += str(breakpoint_info[1]) + '\t'
        row += str(breakpoint_info[2]) + '\t'
        row += str(cycle_info[0]) + '\t'
        row += str(len(cycle_info[1])) + '\t'
        row += (', '.join(cycle_info[1]), 'NA')[len(cycle_info[1]) == 0] + '\t'
        row += (', '.join(dgv_info[cluster_id]), 'NA')[len(dgv_info[cluster_id]) == 0] + '\t'
        row += create_approximate_sequence(sequences, cluster_info[0]) + '\t'
        row += create_exact_sequence(sequences, cluster_info[0], breakpoint_info)
        row += '\n'
        return row


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


    def reverse_complement(sequence):
        return sequence[::-1].translate(string.maketrans('ACTGactg','TGACtgac'))


    def create_approximate_sequence(reference_sequences, cluster_info):
        approximate_sequences = ['', '']
        expected_strands = ('+', '-')
        for side in (0, 1):
            chromosome = cluster_info[side].chromosome
            start = cluster_info[side].start
            end = cluster_info[side].end
            approximate_sequences[side] = reference_sequences[chromosome][start-1:end]
            if cluster_info[side].strand != expected_strands[side]:
                approximate_sequences[side] = reverse_complement(approximate_sequences[side])
        return approximate_sequences[0] + '[]' + approximate_sequences[1]


    def create_exact_sequence(reference_sequences, cluster_info, breakpoint_info):
        if breakpoint_info[0] == 0:
            return 'NA'
        approximate_sequences = ['', '']
        expected_strands = ('+', '-')
        inserted = ''
        if breakpoint_info[0] != 0:
            inserted = 'N' * breakpoint_info[1]
        for side in (0, 1):
            chromosome = cluster_info[side].chromosome
            start = cluster_info[side].start
            end = cluster_info[side].end
            if cluster_info[side].strand == '+':
                end = breakpoint_info[3][side]
            else:
                start = breakpoint_info[3][side]
            approximate_sequences[side] = reference_sequences[chromosome][start-1:end]
            if cluster_info[side].strand != expected_strands[side]:
                approximate_sequences[side] = reverse_complement(approximate_sequences[side])
        return approximate_sequences[0] + '[' + inserted + ']' + approximate_sequences[1]


    def rearrangement_type(cluster):
        if cluster[0].chromosome != cluster[1].chromosome:
            return 'translocation'
        if cluster[0].strand == cluster[1].strand:
            return 'inversion'
        for side1, side2 in ((0, 1), (1, 0)):
            if cluster[side1].strand == '+' and cluster[side2].strand == '-' and cluster[side1].end < cluster[side2].start:
                return 'deletion'
        return 'eversion'


    def filter_clusters(clusters_out_filename, clusters_filename, clusters_prob_filename, align_prob_threshold, chimeric_prob_threshold, valid_prob_threshold, coverage_threshold, readcount_threshold):
        filtered_clusters = set()
        with open(clusters_prob_filename, 'r') as clusters_prob_file:
            for row in csv.reader(clusters_prob_file, delimiter='\t'):
                cluster_id = row[0]
                align_prob = float(row[1])
                chimeric_prob = float(row[2])
                valid_prob = float(row[3])
                if align_prob < align_prob_threshold:
                    continue
                if chimeric_prob < chimeric_prob_threshold:
                    continue
                if valid_prob < valid_prob_threshold:
                    continue
                filtered_clusters.add(cluster_id)
        with open(clusters_filename, 'r') as clusters_file, open(clusters_out_filename, 'w') as clusters_out_file:
            for cluster_id, regions, lib_counts, rows in read_cluster_regions(clusters_file):
                if cluster_id not in filtered_clusters:
                    continue
                if any((regions[cluster_end].end - regions[cluster_end].start + 1 < coverage_threshold for cluster_end in (0, 1))):
                    continue
                if sum(lib_counts.values()) < readcount_threshold:
                    continue
                for row in rows:
                    clusters_out_file.write('\t'.join(row) + '\n')


    class DGVDatabase(object):
        def __init__(self, dgv_filename):
            self.variations = list()
            chrvars = defaultdict(list)
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
        def query(self, regions):
            if regions[0].chromosome != regions[1].chromosome:
                return
            if regions[0].chromosome not in self.intervals:
                return
            break1 = (regions[0].start, regions[0].end)[regions[0].strand == '+']
            break2 = (regions[1].start, regions[1].end)[regions[1].strand == '+']
            start, end = sorted([break1, break2])
            idxs = self.intervals[regions[0].chromosome].find_overlapping(start, end)
            for idx in [idxs[a] for a in range(0, len(idxs))]:
                startdiff = abs(start - self.variations[idx][1])
                enddiff = abs(end - self.variations[idx][2])
                if startdiff < 500 and enddiff < 500:
                    yield self.variations[idx][0]


    def multilib_tabulate(clusters_table, clusters_filename, clusters_prob_filename, clusters_all_filename, breakpoints_filename, genome_fasta, gtf_filename, dgv_filename, cycles_filename, stats):
        cluster_infos = dict()
        lib_names = set(stats.keys())
        with open(clusters_filename, 'r') as clusters_file:
            for cluster_id, regions, lib_counts, rows in read_cluster_regions(clusters_file):
                lib_names.update(lib_counts.keys())
                cluster_infos[cluster_id] = (regions, lib_counts)
        probs = dict()
        with open(clusters_prob_filename, 'r') as clusters_prob_file:
            for row in csv.reader(clusters_prob_file, delimiter='\t'):
                cluster_id = row[0]
                align_prob = float(row[1])
                chimeric_prob = float(row[2])
                valid_prob = float(row[3])
                probs[row[0]] = (align_prob, chimeric_prob, valid_prob)
        breakpoint_infos = dict()
        with open(breakpoints_filename, 'r') as breakpoints_file:
            for cluster_id, breakpoints, num_inserted, split_count, homology in read_breakpoints(breakpoints_file):
                if cluster_id not in cluster_infos:
                    continue
                breakpoint_infos[cluster_id] = (split_count, num_inserted, homology, breakpoints)
        dgv = DGVDatabase(dgv_filename)
        dgv_info = dict()
        for cluster_id, cluster_info in cluster_infos.iteritems():
            dgv_info[cluster_id] = list(dgv.query(cluster_info[0]))
        setcover_ratios = dict()
        with open(clusters_all_filename, 'r') as clusters_all_file:
            for cluster_id, regions, lib_counts, rows in read_cluster_regions(clusters_all_file):
                if cluster_id in cluster_infos:
                    setcover_ratios[cluster_id] = float(sum(cluster_infos[cluster_id][1].values())) / float(sum(lib_counts.values()))
        cycle_infos = dict()
        with open(cycles_filename, 'r') as cycles_file:
            for row in csv.reader(cycles_file, delimiter='\t'):
                score = row[0]
                ids = row[1::4]
                cycle_infos[ids[0]] = (score, ids)
        sequences = dict()
        for id, seq in read_sequences(open(genome_fasta, 'r')):
            sequences[id] = seq
        gene_models = pygenes.GeneModels()
        gene_models.load_ensembl_gtf(gtf_filename)
        with open(clusters_table, 'w') as clusters_table_file:
            clusters_table_file.write(clusters_table_header(lib_names))
            for cluster_id, cluster_info in cluster_infos.iteritems():
                setcover_ratio = setcover_ratios[cluster_id]
                breakpoint_info = breakpoint_infos.get(cluster_id, (0, 0, 0, None))
                cycle_info = cycle_infos.get(cluster_id, ('NA', []))
                clusters_table_file.write(clusters_table_row(lib_names, cluster_id, cluster_info, probs, setcover_ratio, breakpoint_info, dgv_info, cycle_info, sequences, gene_models))


    def merge_tars(output_filename, *input_filename_sets):
        with tarfile.open(output_filename, 'w') as output_tar:
            for input_filenames in input_filename_sets:
                for input_filename in input_filenames.itervalues():
                    with tarfile.open(input_filename, 'r') as in_tar:
                        for tarinfo in in_tar:
                            output_tar.addfile(tarinfo, in_tar.extractfile(tarinfo))


