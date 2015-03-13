import os
import collections
import bisect
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.metrics import roc_curve, auc
import seaborn

import pygenes
import pypeliner

import wrappers
import utils.download


destruct_directory = os.path.abspath(os.path.dirname(__file__))
tools_directory = os.path.join(destruct_directory, 'tools')


def read_simulation_params(sim_config_filename):

    sim_info = {}
    execfile(sim_config_filename, {}, sim_info)

    return sim_info


class BreakpointDatabase(object):
    def __init__(self, breakpoints):
        self.positions = collections.defaultdict(list)
        self.break_ids = collections.defaultdict(set)
        for idx, row in breakpoints.iterrows():
            for side in ('1', '2'):
                self.positions[(row['chromosome_'+side], row['strand_'+side])].append(row['position_'+side])
                self.break_ids[(row['chromosome_'+side], row['strand_'+side], row['position_'+side])].add((row['break_id'], side))
        for key in self.positions.iterkeys():
            self.positions[key] = sorted(self.positions[key])
    def query(self, row, extend=0):
        matched_ids = list()
        for side in ('1', '2'):
            chrom_strand_positions = self.positions[(row['chromosome_'+side], row['strand_'+side])]
            idx = bisect.bisect_left(chrom_strand_positions, row['position_'+side] - extend)
            side_matched_ids = list()
            while idx < len(chrom_strand_positions):
                pos = chrom_strand_positions[idx]
                dist = abs(pos - row['position_'+side])
                if pos >= row['position_'+side] - extend and pos <= row['position_'+side] + extend:
                    for break_id in self.break_ids[(row['chromosome_'+side], row['strand_'+side], pos)]:
                        side_matched_ids.append((break_id, dist))
                if pos > row['position_'+side] + extend:
                    break
                idx += 1
            matched_ids.append(side_matched_ids)
        matched_ids_bypos = list()
        for matched_id_1, dist_1 in matched_ids[0]:
            for matched_id_2, dist_2 in matched_ids[1]:
                if matched_id_1[0] == matched_id_2[0] and matched_id_1[1] != matched_id_2[1]:
                    matched_ids_bypos.append((dist_1 + dist_2, matched_id_1[0]))
        if len(matched_ids_bypos) == 0:
            return None
        return sorted(matched_ids_bypos)[0][1]


def create_tool_wrappers(install_directory):

    tool_wrappers = dict()

    for tool_name, ToolWrapper in wrappers.catalog.iteritems():

        tool_wrappers[tool_name] = ToolWrapper(os.path.join(install_directory, tool_name))

    return tool_wrappers


def run_tool(tool_wrapper, temp_directory, results_filename, control_id=None, **bam_filenames):

    tool_wrapper.run(bam_filenames, results_filename, temp_directory, control_id=control_id)
    

def create_roc_plot(sim_info, tool_wrapper, simulated_filename, predicted_filename, annotated_filename, identified_filename, plot_filename):

    with PdfPages(plot_filename) as pdf:

        breakpoints = pd.read_csv(simulated_filename, sep='\t', header=None,
                          converters={'break_id':str, 'chromosome_1':str, 'chromosome_2':str},
                          names=['break_id',
                                 'chromosome_1', 'strand_1', 'position_1',
                                 'chromosome_2', 'strand_2', 'position_2',
                                 'inserted', 'homology'])

        breakpoints_db = BreakpointDatabase(breakpoints)

        results = pd.read_csv(predicted_filename, sep='\t',
                              converters={'prediction_id':str, 'chromosome_1':str, 'chromosome_2':str})

        min_dist = 200

        def identify_true_positive(row):
            return breakpoints_db.query(row, min_dist)

        results['true_pos_id'] = results.apply(identify_true_positive, axis=1)

        results.to_csv(annotated_filename, sep='\t', index=False, na_rep='NA')

        features = tool_wrapper.features
        invalid_value = -1e9

        identified = breakpoints[['break_id']]
        identified = identified.merge(results[['prediction_id', 'true_pos_id'] + features], left_on='break_id', right_on='true_pos_id', how='outer')
        identified = identified.drop('true_pos_id', axis=1)

        identified.to_csv(identified_filename, sep='\t', index=False, na_rep='NA')

        num_positive = int(sim_info['num_breakpoints'])

        y_test = results['true_pos_id'].notnull()

        tp_identified = np.sum(y_test)

        num_missed = num_positive - np.sum(y_test)

        y_test = np.concatenate([y_test, np.array([True]*num_missed)])

        if y_test.sum() == len(y_test):
            y_test = np.concatenate([y_test, np.array([False])])

        fig = plt.figure(figsize=(16,16))

        for feature in features:
            
            values = results[feature].values
            values = np.concatenate([values, np.array([invalid_value] * (len(y_test) - len(values)))])
            
            tp_counts = list()
            fp_counts = list()
            
            for threshold in np.sort(np.unique(values)):
                
                tp_counts.append(np.sum(y_test[values > threshold]))
                fp_counts.append(np.sum(1 - y_test[values > threshold]))
            
            plt.plot(fp_counts, tp_counts, label=feature.replace('_', ' '))

        plt.title('Features for {0} identified true positives'.format(tp_identified))
        plt.xlabel('False Positive Count')
        plt.ylabel('True Positive Count')
        plt.legend(loc="lower right")

        pdf.savefig(fig)

        for feature in features:

            fig = plt.figure(figsize=(16,16))

            true_features = results.loc[results['true_pos_id'].notnull(), feature].values.astype(float)
            false_features = results.loc[results['true_pos_id'].isnull(), feature].values.astype(float)

            def add_optional_noise(feature_values):
                if len(feature_values) == 1:
                    feature_values = np.concatenate([feature_values, feature_values])
                if len(np.unique(feature_values)) <= 1:
                    feature_values += np.random.randn(*feature_values.shape) * 1e-4
                return feature_values

            true_features = add_optional_noise(true_features)
            false_features = add_optional_noise(false_features)

            seaborn.violinplot([true_features, false_features],
                               names=['True', 'False'])

            plt.title('True vs false for ' + feature)

            pdf.savefig(fig)


def create_genome(chromosomes, include_nonchromosomal, genome_fasta):

    utils.download.download_genome_fasta(genome_fasta,
                                         chromosomes,
                                         include_nonchromosomal)

    pypeliner.commandline.execute('bwa', 'index', genome_fasta)
    pypeliner.commandline.execute('samtools', 'faidx', genome_fasta)

    for extension in ('.fai', '.amb', '.ann', '.bwt', '.pac', '.sa'):
        os.rename(genome_fasta+extension, genome_fasta[:-4]+extension)


def partition_bam(original_filename, output_a_filename, output_b_filename, fraction_a):

    pypeliner.commandline.execute(os.path.join(tools_directory, 'bampartition'),
                                  '-i', original_filename,
                                  '-a', output_a_filename,
                                  '-b', output_b_filename,
                                  '-f', fraction_a)

    pypeliner.commandline.execute('samtools', 'index', output_a_filename, output_a_filename[:-4]+'.bai')
    pypeliner.commandline.execute('samtools', 'index', output_b_filename, output_b_filename[:-4]+'.bai')


def samtools_sort_index(input_filename, output_filename):

    pypeliner.commandline.execute('samtools', 'sort', input_filename, output_filename)

    os.rename(output_filename + '.bam', output_filename)

    assert output_filename.endswith('.tmp')
    index_filename = output_filename[:-4] + '.bai'

    pypeliner.commandline.execute('samtools', 'index', output_filename, index_filename)


def samtools_merge_sort_index(output_filename, *input_filenames):

    sorted_filenames = list()

    for input_filename in input_filenames:
        pypeliner.commandline.execute('samtools', 'sort', input_filename, input_filename+'.sorted')
        sorted_filenames.append(input_filename+'.sorted.bam')

    pypeliner.commandline.execute('samtools', 'merge', output_filename, *sorted_filenames)

    for sorted_filename in sorted_filenames:
        os.remove(sorted_filename)

    assert output_filename.endswith('.tmp')
    index_filename = output_filename[:-4] + '.bai'

    pypeliner.commandline.execute('samtools', 'index', output_filename, index_filename)





