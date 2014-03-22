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
from collections import *

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
    sch.commandline('prepsample', axes, medmem, 'cat', sch.ifile('sample1', axes), sch.ifile('sample2', axes), '>', sch.ofile('sample', axes))
    sch.commandline('bwtsamplek2', axes, medmem, cfg.bowtie2_bin, '--very-sensitive', '-k', '2', '-x', cfg.genome_fasta, sch.ifile('sample', axes), '>', sch.ofile('sample.k2.sam', axes))
    sch.commandline('alignnull', axes, medmem, cfg.alignnull_tool, '-g', cfg.gap_score, '-x', cfg.mismatch_score, '-m', cfg.match_score, '-c', sch.ifile('sample.sam', axes), '-a', sch.ifile('sample.k2.sam', axes), '-s', sch.ifile('sample', axes), '-r', cfg.genome_fasta, '>', sch.ofile('samples.align.null', axes))
    sch.transform('scorestats', axes, medmem, score_stats.create_score_stats, None, sch.ifile('samples.align.true', axes), sch.ifile('samples.align.null', axes), int(cfg.match_score), sch.ofile('score.stats', axes))

    sch.transform('readstats', axes, medmem, calculate_concordant_stats, sch.oobj('stats', axes), sch.ifile('sample.sam', axes))

    sch.transform('prepseed', axes, medmem, prepare_seed_fastq, None, sch.ifile('reads1', axes), sch.ifile('reads2', axes), 36, sch.ofile('reads.seed', axes))
    sch.commandline('prepreads', axes, medmem, 'cat', sch.ifile('reads1', axes), sch.ifile('reads2', axes), '>', sch.ofile('reads', axes))
    sch.commandline('bwtseed', axes, medmem, cfg.bowtie_bin, cfg.genome_fasta, sch.ifile('reads.seed', axes), '--chunkmbs', '512', '-k', '1000', '-m', '1000', '--strata', '--best', '-S', '>', sch.ofile('reads.seed.sam', axes))
    sch.commandline('realign', axes, medmem, cfg.realign2_tool, '-a', sch.ifile('reads.seed.sam', axes), '-s', sch.ifile('reads', axes), '-r', cfg.genome_fasta, '-g', cfg.gap_score, '-x', cfg.mismatch_score, '-m', cfg.match_score, '--flmin', sch.iobj('stats', axes).prop('fragment_length_min'), '--flmax', sch.iobj('stats', axes).prop('fragment_length_max'), '--tchimer', cfg.chimeric_threshold, '--talign', cfg.alignment_threshold, '--pchimer', cfg.chimeric_prior, '--pvalid', cfg.readvalid_prior, '--tvalid', cfg.readvalid_threshold, '-z', sch.ifile('score.stats', axes), '--span', sch.ofile('spanning.alignments', axes), '--split', sch.ofile('split.alignments', axes))

    sch.commandline('cluster', axes, himem, cfg.clustermatepairs_tool, '-a', sch.ifile('spanning.alignments', axes), '-m', '1', '-u', sch.iobj('stats', axes).prop('fragment_length_mean'), '-d', sch.iobj('stats', axes).prop('fragment_length_stddev'), '-o', sch.ofile('clusters', axes))
    sch.commandline('breaks', axes, himem, cfg.predictbreaks_tool, '-s', sch.ifile('split.alignments', axes), '-c', sch.ifile('clusters', axes), '-r', cfg.genome_fasta, '-b', sch.ofile('breakpoints', axes))
    sch.transform('results', axes, lowmem, compile_results, None, sch.ifile('simulated.info', axes), sch.ifile('breakpoints', axes), sch.ifile('clusters', axes), sch.ofile('identified', axes), sch.ofile('classify', axes))

    sch.transform('collate', (), lowmem, collate_results, None, sch.iobj('simulation.params', axes), sch.ifile('identified', axes), sch.ifile('classify', axes), sch.output(args.results))

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


    def collate_results(sim_params, identified_filenames, classify_filenames, results_filename):
        with open(results_filename, 'w') as results_file:
            fields = list(set([field for sim in sim_params.values() for field in sim.keys()]))
            results_file.write('\t'.join(['sim_id'] + fields + ['sensitivity_exact', 'sensitivity_approx', 'homology_sensitivity', 'specificity']) + '\n')
            for sim_id in sim_params.keys():
                results_file.write(str(sim_id) + '\t')
                results_file.write('\t'.join((sim_params[sim_id][field] for field in fields)))
                with open(identified_filenames[sim_id], 'r') as identified_file:
                    total = 0.0
                    identified_exact = 0.0
                    identified_approx = 0.0
                    homology_identified = 0.0
                    for row in csv.reader(identified_file, delimiter='\t'):
                        total += 1.0
                        if row[1] != 'N':
                            identified_exact += 1.0
                        if row[2] != 'N':
                            identified_approx += 1.0
                        if row[3] != 'N':
                            homology_identified += 1.0
                    sensitivity_exact = identified_exact / total
                    sensitivity_approx = identified_approx / total
                    homology_sensitivity = homology_identified / identified_exact
                with open(classify_filenames[sim_id], 'r') as classify_file:
                    total_clusters = 0.0
                    total_true = 0.0
                    for row in csv.reader(classify_file, delimiter='\t'):
                        total_clusters += 1.0
                        if row[1] == '1':
                            total_true += 1.0
                    specificity = total_true / total_clusters
                results_file.write('\t' + '\t'.join([str(sensitivity_exact), str(sensitivity_approx), str(homology_sensitivity), str(specificity)]) + '\n')


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
            region = [[None, None, [], []], [None, None, [], []]]
            for row in rows:
                cluster_end = int(row[1])
                region[cluster_end][0] = row[4]
                region[cluster_end][1] = row[5]
                region[cluster_end][2].append(int(row[6]))
                region[cluster_end][3].append(int(row[7]))
            for cluster_end in (0, 1):
                region[cluster_end][2] = min(region[cluster_end][2])
                region[cluster_end][3] = max(region[cluster_end][3])
            yield cluster_id, region


    def match(region, breakpoint):
        min_dist = 200
        for side in (0, 1):
            if region[0][0] != breakpoint[side][0] or region[0][1] != breakpoint[side][1] or region[1][0] != breakpoint[1-side][0] or region[1][1] != breakpoint[1-side][1]:
                continue
            if breakpoint[side][2] > region[0][3] + min_dist or breakpoint[side][2] < region[0][2] - min_dist or breakpoint[1-side][2] > region[1][3] + min_dist or breakpoint[1-side][2] < region[1][2] - min_dist:
                continue
            return True
        return False


    def compile_results(sim_filename, breakpoints_filename, clusters_filename, identified_filename, classify_filename):
        predictions = defaultdict(list)
        with open(breakpoints_filename, 'r') as breakpoints_file:
            for prediction_id, cluster_id, breakends in read_breakends(breakpoints_file):
                predictions[breakends].append(cluster_id)
        predicted_homology = dict()
        with open(breakpoints_filename, 'r') as breakpoints_file:
            for cluster_id, homology in read_homology(breakpoints_file):
                predicted_homology[cluster_id] = homology
        cluster_regions = dict()
        with open(clusters_filename, 'r') as clusters_file:
            for cluster_id, region in read_cluster_regions(clusters_file):
                cluster_regions[cluster_id] = region
        with open(sim_filename, 'r') as sim_file, open(identified_filename, 'w') as identified_file, open(classify_filename, 'w') as classify_file:
            true_cluster_ids = set()
            for sim_id, breakends, homology in read_sim(sim_file):
                exact_cluster_id = 'N'
                approx_cluster_id = 'N'
                homology_correct = 'N'
                for cluster_id in predictions[breakends]:
                    exact_cluster_id = cluster_id
                    true_cluster_ids.add(cluster_id)
                    if predicted_homology[cluster_id] == homology:
                        homology_correct = 'Y'
                for cluster_id, region in cluster_regions.iteritems():
                    if match(region, list(breakends)):
                        approx_cluster_id = cluster_id
                        true_cluster_ids.add(approx_cluster_id)
                identified_file.write('\t'.join((str(sim_id), str(exact_cluster_id), str(approx_cluster_id), homology_correct)) + '\n')
            for cluster_id, region in cluster_regions.iteritems():
                classify_file.write('\t'.join([str(cluster_id), str(int(cluster_id in true_cluster_ids))]) + '\n')


    class ConcordantReadStatsAccumulator(object):
        def __init__(self):
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
            return int(max(0, self.fragment_length_mean - 3 * self.fragment_length_stddev))
        @property
        def fragment_length_max(self):
            return int(self.fragment_length_mean + 3 * self.fragment_length_stddev)


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


    def calculate_concordant_stats(concordant_sam):
        stats = ConcordantReadStatsAccumulator()
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
