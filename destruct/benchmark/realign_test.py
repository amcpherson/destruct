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
import collections
import pandas as pd

import pypeliner

if __name__ == '__main__':

    import realign_test
    from realign_test import *

    import destruct
    from destruct import *

    import create_breakpoint_simulation

    import score_stats

    parser = argparse.ArgumentParser()
    pypeliner.easypypeliner.add_arguments(parser)
    parser.add_argument('simconfig', help='Simulation Configuration filename')
    parser.add_argument('results', help='Results filename')
    args = parser.parse_args()

    cfg = pypeliner.easypypeliner.Config(vars(parser.parse_args()))
    pyp = pypeliner.easypypeliner.EasyPypeliner([destruct, realign_test], cfg)
    sch = pyp.sch

    lowmem = {'mem':1}
    medmem = {'mem':8}
    himem = {'mem':32}

    axes = ('bysim',)

    sch.transform('create', (), lowmem, create_simulation_params, sch.oobj('simulation.params', axes), sch.input(args.simconfig))

    sch.transform('generate', axes, medmem, create_breakpoint_simulation.create, None, cfg, sch.iobj('simulation.params', axes), sch.ofile('simulated.fasta', axes), sch.ofile('simulated.info', axes), sch.ofile('concordant1', axes), sch.ofile('concordant2', axes), sch.ofile('discordant1', axes), sch.ofile('discordant2', axes))

    sch.transform('index1', axes, lowmem, index_fastq, None, sch.ifile('concordant1', axes), sch.ofile('sample1', axes))
    sch.transform('index2', axes, lowmem, index_fastq, None, sch.ifile('concordant2', axes), sch.ofile('sample2', axes))
    sch.transform('index3', axes, lowmem, index_fastq, None, sch.ifile('discordant1', axes), sch.ofile('reads1', axes))
    sch.transform('index4', axes, lowmem, index_fastq, None, sch.ifile('discordant2', axes), sch.ofile('reads2', axes))

    sch.commandline('bwtsample', axes, medmem, cfg.bowtie2_bin, '--no-mixed', '--no-discordant', '--very-sensitive', '-X', cfg.fragment_length_max, cfg.bowtie2_options, '-x', cfg.genome_fasta, '-1', sch.ifile('sample1', axes), '-2', sch.ifile('sample2', axes), '>', sch.ofile('sample.sam', axes))
    sch.commandline('aligntrue', axes, medmem, cfg.aligntrue_tool, '-g', cfg.gap_score, '-x', cfg.mismatch_score, '-m', cfg.match_score, '-r', cfg.genome_fasta, '-a', sch.ifile('sample.sam', axes), '>', sch.ofile('samples.align.true', axes))
    sch.transform('scorestats', axes, medmem, score_stats.create_score_stats, None, sch.ifile('samples.align.true', axes), int(cfg.match_score), sch.ofile('score.stats', axes), sch.ofile('score.stats.plots', axes), sch.inst('bysim'))

    sch.transform('readstats', axes, medmem, calculate_concordant_stats, sch.oobj('stats', axes), sch.ifile('sample.sam', axes), float(cfg.fragment_length_num_stddevs))

    sch.transform('prepseed', axes, medmem, prepare_seed_fastq, None, sch.ifile('reads1', axes), sch.ifile('reads2', axes), 36, sch.ofile('reads.seed', axes))
    sch.commandline('prepreads', axes, medmem, 'cat', sch.ifile('reads1', axes), sch.ifile('reads2', axes), '>', sch.ofile('reads', axes))
    sch.commandline('bwtseed', axes, medmem, cfg.bowtie_bin, cfg.genome_fasta, sch.ifile('reads.seed', axes), '--chunkmbs', '512', '-k', '1000', '-m', '1000', '--strata', '--best', '-S', '>', sch.ofile('reads.seed.sam', axes))
    sch.commandline('realign', axes, medmem, cfg.realign2_tool, '-a', sch.ifile('reads.seed.sam', axes), '-s', sch.ifile('reads', axes), '-r', cfg.genome_fasta, '-g', cfg.gap_score, '-x', cfg.mismatch_score, '-m', cfg.match_score, '--flmin', sch.iobj('stats', axes).prop('fragment_length_min'), '--flmax', sch.iobj('stats', axes).prop('fragment_length_max'), '--tchimer', cfg.chimeric_threshold, '--talign', cfg.alignment_threshold, '--pchimer', cfg.chimeric_prior, '--tvalid', cfg.readvalid_threshold, '-z', sch.ifile('score.stats', axes), '--span', sch.ofile('spanning.alignments', axes), '--split', sch.ofile('split.alignments', axes))

    sch.transform('eval_span_align', axes, medmem, eval_span_align, None, sch.ifile('simulated.info', axes), sch.ifile('reads1', axes), sch.ifile('spanning.alignments', axes), sch.ofile('spanning.alignments.anno', axes), sch.ofile('spanning.alignments.eval', axes))
    sch.transform('eval_split_align', axes, medmem, eval_split_align, None, sch.ifile('simulated.info', axes), sch.ifile('reads1', axes), sch.ifile('split.alignments', axes), sch.ofile('split.alignments.anno', axes), sch.ofile('split.alignments.eval', axes))

    sch.commandline('cluster', axes, himem, cfg.clustermatepairs_tool, '-a', sch.ifile('spanning.alignments', axes), '-m', '1', '-u', sch.iobj('stats', axes).prop('fragment_length_mean'), '-d', sch.iobj('stats', axes).prop('fragment_length_stddev'), '-o', sch.ofile('clusters', axes))
    sch.commandline('breaks', axes, himem, cfg.predictbreaks_tool, '-s', sch.ifile('split.alignments', axes), '-c', sch.ifile('clusters', axes), '-r', cfg.genome_fasta, '-b', sch.ofile('breakpoints', axes))

    sch.transform('results', axes, lowmem, compile_results, None, sch.ifile('simulated.info', axes), sch.ifile('breakpoints', axes), sch.ifile('clusters', axes), sch.ifile('spanning.alignments.eval', axes), sch.ifile('split.alignments.eval', axes), sch.ofile('identified', axes), sch.ofile('classified', axes))

    sch.transform('collate', (), lowmem, collate_results, None, sch.iobj('simulation.params', axes), sch.ifile('identified', axes), sch.ifile('classified', axes), sch.output(args.results))

    pyp.run()

else:


    def create_simulation_params(config_filename):
        config = ConfigParser.ConfigParser()
        config.read(config_filename)
        section = 'main'
        sim_info_base = dict()
        sim_info_base['adjacent_length'] = config.get(section, 'adjacent_length')
        sim_info_base['breakpoints_seed'] = config.get(section, 'breakpoints_seed')
        sim_info_base['dwgsim_seed'] = config.get(section, 'dwgsim_seed')
        sim_info_base['num_breakpoints'] = config.get(section, 'num_breakpoints')
        sim_info_base['num_concordant'] = config.get(section, 'num_concordant')
        read_length_list = config.get(section, 'read_length_list').split()
        fragment_mean_list = config.get(section, 'fragment_mean_list').split()
        fragment_stddev_list = config.get(section, 'fragment_stddev_list').split()
        coverage_list = config.get(section, 'coverage_list').split()
        num_inserted_list = config.get(section, 'num_inserted_list').split()
        homology_list = config.get(section, 'homology_list').split()
        sim_id = 0
        sim_params = dict()
        for read_length, fragment_mean, fragment_stddev, coverage, num_inserted, homology in itertools.product(read_length_list, fragment_mean_list, fragment_stddev_list, coverage_list, num_inserted_list, homology_list):
            sim_params[sim_id] = dict(sim_info_base)
            sim_params[sim_id]['read_length'] = read_length
            sim_params[sim_id]['fragment_mean'] = fragment_mean
            sim_params[sim_id]['fragment_stddev'] = fragment_stddev
            sim_params[sim_id]['coverage'] = coverage
            sim_params[sim_id]['num_inserted'] = num_inserted
            sim_params[sim_id]['homology'] = homology
            sim_id += 1
        return sim_params


    def create_sim_read_table(reads_filename):
        sim_reads = list()
        with open(reads_filename, 'r') as reads_file:
                for name, seq, comment, qual in itertools.izip_longest(*[reads_file]*4):
                        assert name[0] == '@'
                        read_id = int(name[1:].split('/')[0])
                        assert comment[:2] == '+@'
                        sim_id = int(comment[2:].split('_')[0])
                        sim_reads.append((sim_id, read_id))
        sim_reads = pd.DataFrame(sim_reads, columns=['sim_id', 'read_id'])
        return sim_reads


    def read_simulated_table(simulated_filename):
        simulated = pd.read_csv(simulated_filename, sep='\t', header=None,
                                converters={'chromosome_1':str, 'chromosome_2':str},
                                names=['sim_id',
                                       'chromosome_1', 'strand_1', 'position_1',
                                       'chromosome_2', 'strand_2', 'position_2',
                                       'inserted', 'homology'])
        return simulated


    def eval_span_align(simulated_filename, reads_filename, spanning_filename, spanning_annotated_filename, eval_filename):

        max_diff = 1000

        sim_reads = create_sim_read_table(reads_filename)
        simulated = read_simulated_table(simulated_filename)

        spanning = pd.read_csv(spanning_filename, sep='\t', header=None, index_col=False,
                               converters={'chromosome':str},
                               names=['read_id', 'read_end',
                                      'chromosome', 'strand', 'start', 'end',
                                      'read_length', 'matched_length', 'score',
                                      'p_align', 'p_chrimeric', 'p_valid'])

        assert len(spanning.index) > 0

        spanning['approx_break'] = (spanning['strand'] == '+') * spanning['end'] + (spanning['strand'] != '+') * spanning['start']

        spanning.set_index(['read_id', 'read_end'], inplace=True)
        spanning['max_p_align'] = spanning.groupby(level=[0,1])['p_align'].max()
        spanning.reset_index(inplace=True)
        spanning['is_best'] = spanning['max_p_align'] == spanning['p_align']

        spanning = spanning.merge(sim_reads, left_on='read_id', right_on='read_id')
        spanning = spanning.merge(simulated, left_on='sim_id', right_on='sim_id')

        for side in ('1', '2'):
            spanning['is_side_'+side] = (spanning['chromosome'] == spanning['chromosome_'+side]) & \
                                        (spanning['strand'] == spanning['strand_'+side]) & \
                                        (spanning['approx_break'] - max_diff <= spanning['position_'+side]) & \
                                        (spanning['approx_break'] + max_diff >= spanning['position_'+side])

        spanning['is_true'] = spanning['is_side_1'] | spanning['is_side_2']
        spanning['is_true_best'] = spanning['is_true'] & spanning['is_best']

        spanning.to_csv(spanning_annotated_filename, sep='\t', index=False)

        spanning = spanning.set_index(['sim_id', 'read_id', 'read_end'])

        true_alignment = spanning['is_true'].groupby(level=[0,1,2]).apply(any).unstack().apply(any, axis=1)
        true_alignment = true_alignment * 1
        true_alignment = true_alignment.groupby(level=0).mean()

        true_and_best_alignment = spanning['is_true_best'].groupby(level=[0,1,2]).apply(any).unstack().apply(any, axis=1)
        true_and_best_alignment = true_and_best_alignment * 1
        true_and_best_alignment = true_and_best_alignment.groupby(level=0).mean()

        span_eval = pd.concat([true_alignment, true_and_best_alignment], axis=1, keys=['true_align', 'true_best_align'])

        span_eval.to_csv(eval_filename, sep='\t', index=True)


    def eval_split_align(simulated_filename, reads_filename, split_filename, split_annotated_filename, eval_filename):

        sim_reads = create_sim_read_table(reads_filename)
        simulated = read_simulated_table(simulated_filename)

        split = pd.read_csv(split_filename, sep='\t', header=None, index_col=False,
                            converters={'chromosome_1':str, 'chromosome_2':str},
                            names=['read_id', 'read_end',
                                   'chromosome_1', 'strand_1', 'position_1',
                                   'chromosome_2', 'strand_2', 'position_2',
                                   'inserted', 'length_1', 'length_2',
                                   'score_1', 'score_2', 'score', 'p_align'])

        split.set_index(['read_id', 'read_end'], inplace=True)
        split['max_p_align'] = split.groupby(level=[0,1])['p_align'].max()
        split.reset_index(inplace=True)
        split['is_best'] = split['max_p_align'] == split['p_align']

        split = split.merge(sim_reads, left_on='read_id', right_on='read_id')

        merge_columns = ['sim_id']
        for side in ('1', '2'):
            merge_columns.append('chromosome_'+side)
            merge_columns.append('strand_'+side)
            merge_columns.append('position_'+side)

        sim_side_1 = simulated[merge_columns]
        sim_side_1['is_side_1'] = True

        split = split.merge(sim_side_1, left_on=merge_columns, right_on=merge_columns, how='left')

        def flip_side(col):
            if col.endswith('_2'):
                return col[:-2]+'_1'
            elif col.endswith('_1'):
                return col[:-2]+'_2'
            else:
                return col

        sim_side_2 = simulated[merge_columns]
        sim_side_2 = sim_side_2.rename(columns=flip_side)
        sim_side_2['is_side_2'] = True

        split = split.merge(sim_side_2, left_on=merge_columns, right_on=merge_columns, how='left')

        split['is_side_1'] = split['is_side_1'].fillna(False)
        split['is_side_2'] = split['is_side_2'].fillna(False)

        split['is_true'] = split['is_side_1'] | split['is_side_2']
        split['is_true_best'] = split['is_true'] & split['is_best']

        split.to_csv(split_annotated_filename, sep='\t', index=False)

        split = split.set_index(['sim_id', 'read_id', 'read_end'])

        true_alignment = split['is_true'].groupby(level=[0,1,2]).apply(any)
        true_alignment = true_alignment * 1
        true_alignment = true_alignment.groupby(level=0).mean()

        true_and_best_alignment = split['is_true_best'].groupby(level=[0,1,2]).apply(any)
        true_and_best_alignment = true_and_best_alignment * 1
        true_and_best_alignment = true_and_best_alignment.groupby(level=0).mean()

        split_eval = pd.concat([true_alignment, true_and_best_alignment], axis=1, keys=['true_align', 'true_best_align'])

        split_eval.to_csv(eval_filename, sep='\t', index=True)


    def collate_results(sim_params, identified_filenames, classified_filenames, results_filename):

        with open(results_filename, 'w') as results_file:

            fields = list(set([field for sim in sim_params.values() for field in sim.keys()]))
            
            results_file.write('\t'.join(['sim_id'] + fields + ['sensitivity_exact', 'sensitivity_approx', 'homology_sensitivity', 'specificity', 'read_count', 'span_true_align', 'span_true_best_align', 'split_true_align', 'split_true_best_align']) + '\n')
            
            for sim_id in sim_params.keys():

                results_file.write(str(sim_id) + '\t')
                results_file.write('\t'.join((sim_params[sim_id][field] for field in fields)) + '\t')

                identified = pd.read_csv(identified_filenames[sim_id], sep='\t')
                sensitivity_exact = identified['exact_cluster_id'].notnull().mean()
                sensitivity_approx = identified['approx_cluster_id'].notnull().mean()
                homology_sensitivity = identified['homology_correct'].mean()
                read_count = identified['read_count'].mean()
                span_true_align = identified['span_true_align'].mean()
                span_true_best_align = identified['span_true_best_align'].mean()
                split_true_align = identified['split_true_align'].mean()
                split_true_best_align = identified['split_true_best_align'].mean()

                classified = pd.read_csv(classified_filenames[sim_id], sep='\t')
                specificity = classified['is_true'].mean()

                results_file.write(str(sensitivity_exact) + '\t')
                results_file.write(str(sensitivity_approx) + '\t')
                results_file.write(str(homology_sensitivity) + '\t')
                results_file.write(str(specificity) + '\t')
                results_file.write(str(read_count) + '\t')
                results_file.write(str(span_true_align) + '\t')
                results_file.write(str(span_true_best_align) + '\t')
                results_file.write(str(split_true_align) + '\t')
                results_file.write(str(split_true_best_align) + '\n')


    def read_sim(input_file):
        for row in csv.reader(input_file, delimiter='\t'):
            sim_id = int(row[0])
            breakends = frozenset([(row[1],row[2],int(row[3])),(row[4],row[5],int(row[6]))])
            homology = int(row[8])
            yield sim_id, breakends, homology


    def read_breakends(input_file):
        prediction_id = 0
        for row in csv.reader(input_file, delimiter='\t'):
            cluster_id = int(row[0])
            breakends = frozenset([(row[1],row[2],int(row[3])),(row[4],row[5],int(row[6]))])
            yield prediction_id, cluster_id, breakends
            prediction_id += 1


    def read_homology(input_file):
        for cluster_id, rows in itertools.groupby(csv.reader(input_file, delimiter='\t'), lambda row: int(row[0])):
            cluster_id = int(cluster_id)
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
            yield cluster_id, homology


    def read_cluster_regions(input_file):
        for cluster_id, rows in itertools.groupby(csv.reader(input_file, delimiter='\t'), lambda row: row[0]):
            cluster_id = int(cluster_id)
            region = [[None, None, [], []], [None, None, [], []]]
            read_ids = set()
            for row in rows:
                cluster_end = int(row[1])
                region[cluster_end][0] = row[4]
                region[cluster_end][1] = row[5]
                region[cluster_end][2].append(int(row[6]))
                region[cluster_end][3].append(int(row[7]))
                read_ids.add(row[2])
            for cluster_end in (0, 1):
                region[cluster_end][2] = min(region[cluster_end][2])
                region[cluster_end][3] = max(region[cluster_end][3])
            yield cluster_id, region, len(read_ids)


    def match(region, breakpoint):
        min_dist = 200
        for side in (0, 1):
            if region[0][0] != breakpoint[side][0] or region[0][1] != breakpoint[side][1] or region[1][0] != breakpoint[1-side][0] or region[1][1] != breakpoint[1-side][1]:
                continue
            if breakpoint[side][2] > region[0][3] + min_dist or breakpoint[side][2] < region[0][2] - min_dist or breakpoint[1-side][2] > region[1][3] + min_dist or breakpoint[1-side][2] < region[1][2] - min_dist:
                continue
            return True
        return False


    def compile_results(sim_filename, breakpoints_filename, clusters_filename, span_eval_filename, split_eval_filename, identified_filename, classified_filename):

        predictions = collections.defaultdict(list)
        with open(breakpoints_filename, 'r') as breakpoints_file:
            for prediction_id, cluster_id, breakends in read_breakends(breakpoints_file):
                predictions[breakends].append(cluster_id)
        
        predicted_homology = dict()
        with open(breakpoints_filename, 'r') as breakpoints_file:
            for cluster_id, homology in read_homology(breakpoints_file):
                predicted_homology[cluster_id] = homology
        
        cluster_regions = dict()
        cluster_read_counts = dict()
        with open(clusters_filename, 'r') as clusters_file:
            for cluster_id, region, read_count in read_cluster_regions(clusters_file):
                cluster_regions[cluster_id] = region
                cluster_read_counts[cluster_id] = read_count
        
        span_eval = pd.read_csv(span_eval_filename, sep='\t', index_col='sim_id')
        
        split_eval = pd.read_csv(split_eval_filename, sep='\t', index_col='sim_id')
        
        with open(sim_filename, 'r') as sim_file, open(identified_filename, 'w') as identified_file, open(classified_filename, 'w') as classified_file:

            true_cluster_ids = set()
            identified_file.write('\t'.join(['sim_id', 'exact_cluster_id', 'approx_cluster_id', 'homology_correct', 'read_count', 'span_true_align', 'span_true_best_align', 'split_true_align', 'split_true_best_align']) + '\n')

            for sim_id, breakends, homology in read_sim(sim_file):

                exact_cluster_id = 'NA'
                approx_cluster_id = 'NA'
                homology_correct = '0'
                for cluster_id in predictions[breakends]:
                    exact_cluster_id = cluster_id
                    true_cluster_ids.add(cluster_id)
                    if predicted_homology[cluster_id] == homology:
                        homology_correct = '1'
                for cluster_id, region in cluster_regions.items():
                    if match(region, list(breakends)):
                        approx_cluster_id = cluster_id
                        true_cluster_ids.add(approx_cluster_id)

                read_count = 'NA'
                if exact_cluster_id != 'NA':
                    read_count = cluster_read_counts[exact_cluster_id]
                elif approx_cluster_id != 'NA':
                    read_count = cluster_read_counts[approx_cluster_id]

                identified_file.write(str(sim_id) + '\t')
                identified_file.write(str(exact_cluster_id) + '\t')
                identified_file.write(str(approx_cluster_id) + '\t')
                identified_file.write(str(homology_correct) + '\t')
                identified_file.write(str(read_count) + '\t')
                identified_file.write(str(span_eval['true_align'].get(sim_id, default=0.0)) + '\t')
                identified_file.write(str(span_eval['true_best_align'].get(sim_id, default=0.0)) + '\t')
                identified_file.write(str(split_eval['true_align'].get(sim_id, default=0.0)) + '\t')
                identified_file.write(str(split_eval['true_best_align'].get(sim_id, default=0.0)) + '\n')
            
            classified_file.write('\t'.join(['cluster_id', 'is_true']) + '\n')
            for cluster_id, region in cluster_regions.items():
                classified_file.write('\t'.join([str(cluster_id), str(int(cluster_id in true_cluster_ids))]) + '\n')


    class ConcordantReadStatsAccumulator(object):
        def __init__(self, fragment_length_num_stddevs):
            self.fragment_length_num_stddevs = fragment_length_num_stddevs
            self.read_lengths = set()
            self.fragment_length_sum = 0
            self.fragment_length_sum_sq = 0
            self.fragment_count = 0
        def add(self, read_length_1, read_length_2, fragment_length):
            self.read_lengths.add(read_length_1)
            self.read_lengths.add(read_length_2)
            self.fragment_length_sum += float(fragment_length)
            self.fragment_length_sum_sq += float(fragment_length) ** 2
            self.fragment_count += 1
        def update(self, other):
            self.read_lengths.update(other.read_lengths)
            self.fragment_length_sum += other.fragment_length_sum
            self.fragment_length_sum_sq += other.fragment_length_sum_sq
            self.fragment_count += other.fragment_count
        @property
        def fragment_length_mean(self):
            return self.fragment_length_sum / self.fragment_count
        @property
        def fragment_length_stddev(self):
            return self.fragment_length_variance ** 0.5
        @property
        def fragment_length_variance(self):
            return self.fragment_length_sum_sq / self.fragment_count - self.fragment_length_mean ** 2
        @property
        def fragment_length_min(self):
            return int(max(0, self.fragment_length_mean - self.fragment_length_num_stddevs * self.fragment_length_stddev))
        @property
        def fragment_length_max(self):
            return int(self.fragment_length_mean + self.fragment_length_num_stddevs * self.fragment_length_stddev)


    class SamRecord(object):
        def __init__(self, line):
            self.line = line
            self.fields = line.rstrip().split('\t')
            self.fields[1] = int(self.fields[1])
        @property
        def is_read_mapped(self):
            return self.fields[1] & 0x04 == 0
        @property
        def is_pair_mapped(self):
            return self.fields[1] & 0x02 != 0
        @property
        def read_id(self):
            return self.fields[0]
        @property
        def fragment_id(self):
            return self.fields[0]
        @property
        def read_end(self):
            if self.fields[1] & 0x40:
                return 0
            elif self.fields[1] & 0x80:
                return 1
            else:
                return None
        @property
        def read_length(self):
            return len(self.fields[9])
        @property
        def chromosome(self):
            return self.fields[2]
        @property
        def position(self):
            return int(self.fields[3])
        @property
        def strand(self):
            if self.fields[1] & 0x10:
                return '-'
            else:
                return '+'
        @property
        def isize(self):
            return abs(int(self.fields[8]))
        @property
        def alignment_score(self):
            for opt in self.fields[11:]:
                if opt.startswith('AS:i:'):
                    return int(opt[5:])
        @property
        def num_matches(self):
            return self.matched - self.num_mismatches
        @property
        def num_mismatches(self):
            for opt in self.fields[11:]:
                if opt.startswith('XM:i:'):
                    return int(opt[5:])
        @property
        def start_gap(self):
            start_gap_match = self.start_gap_re.search(self.fields[5])
            if start_gap_match is not None:
                return int(start_gap_match.group(1))
            else:
                return 0
        start_gap_re = re.compile('^(\d+)S')
        @property
        def end_gap(self):
            end_gap_match = self.end_gap_re.search(self.fields[5])
            if end_gap_match is not None:
                return int(end_gap_match.group(1))
            else:
                return 0
        end_gap_re = re.compile('(\d+)S$')
        @property
        def matched(self):
            matched = 0
            for count in self.matched_re.findall(self.fields[5]):
                matched += int(count)
            return matched
        matched_re = re.compile('(\d+)M')
        @property
        def deleted(self):
            deleted = 0
            for count in self.deleted_re.findall(self.fields[5]):
                deleted += int(count)
            return deleted
        deleted_re = re.compile('(\d+)D')


    def calculate_concordant_stats(concordant_sam, fragment_length_num_stddevs):
        stats = ConcordantReadStatsAccumulator(fragment_length_num_stddevs)
        try:
            with open(concordant_sam, 'r') as concordant:
                for line_1 in concordant:
                    if line_1[0] == '@':
                        continue
                    line_2 = concordant.next()
                    record_1 = SamRecord(line_1)
                    record_2 = SamRecord(line_2)
                    if (record_1.fragment_id != record_2.fragment_id):
                        logging.getLogger('pipeline').error('Mismatched read names %s and %s' % (record_1.read_id, record_2.read_id))
                    if record_1.isize != record_2.isize:
                        logging.getLogger('pipeline').error('Mismatched isizes %d and %d for reads %s and %s' % (record_1.isize, record_2.isize, record_1.read_id, record_2.read_id))
                    if record_1.chromosome != record_2.chromosome:
                        logging.getLogger('pipeline').error('Mismatched references for concordant reads %s and %s' % (record_1.read_id, record_2.read_id))
                    if not record_1.is_pair_mapped or not record_2.is_pair_mapped:
                        continue
                    stats.add(record_1.read_length, record_2.read_length, record_1.isize)
        except StopIteration:
            pass
        return stats


    def index_fastq(in_filename, out_filename):
        with open(in_filename, 'r') as in_file, open(out_filename, 'w') as out_file:
            read_index = 0
            for name, seq, comment, qual in itertools.izip_longest(*[in_file]*4):
                out_file.write('@' + str(read_index) + name.rstrip()[-2:] + '\n')
                out_file.write(seq)
                out_file.write('+' + name)
                out_file.write(qual)
                read_index += 1
