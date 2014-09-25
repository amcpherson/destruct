import csv
import sys
import logging
import os
import ConfigParser
import re
import itertools
import collections
import bisect
import subprocess
import argparse
import string
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.metrics import roc_curve, auc
import seaborn

import pygenes
import pypeliner
import pypeliner.managed as mgd

import wrappers
import utils.download

destruct_directory = os.path.abspath(os.path.dirname(__file__))

tools_directory = os.path.join(destruct_directory, 'tools')
default_config_filename = os.path.join(destruct_directory, 'defaultconfig.py')


if __name__ == '__main__':

    import destruct_bam_test
    import create_breakpoint_simulation

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)
    argparser.add_argument('simconfig', help='Simulation configuration filename')
    argparser.add_argument('bam', help='Source bam filename')
    argparser.add_argument('ref', help='Reference genome for source bam')
    argparser.add_argument('installdir', help='Tool installations directory')
    argparser.add_argument('outdir', help='Output directory')
    argparser.add_argument('-c', '--config', help='Configuration filename')

    args = vars(argparser.parse_args())

    config = {}

    if args['config'] is not None:
        execfile(args['config'], {}, config)

    config.update(args)

    pyp = pypeliner.app.Pypeline([destruct_bam_test, create_breakpoint_simulation], config)

    try:
        os.makedirs(args['outdir'])
    except OSError:
        pass

    ctx = {'mem':4}

    pyp.sch.transform('read_params', (), ctx,
        destruct_bam_test.read_simulation_params,
        mgd.TempOutputObj('simulation.params'),
        mgd.InputFile(args['simconfig']))

    pyp.sch.transform('create_sim', (), ctx,
        create_breakpoint_simulation.create_breakpoints,
        None,
        mgd.TempInputObj('simulation.params'),
        mgd.InputFile(args['ref']),
        mgd.OutputFile(os.path.join(args['outdir'], 'simulated.fasta')),
        mgd.OutputFile(os.path.join(args['outdir'], 'simulated.tsv')))

    pyp.sch.transform('partition', (), ctx,
        destruct_bam_test.partition_bam,
        None,
        mgd.InputFile(args['bam']),
        mgd.OutputFile(os.path.join(args['outdir'], 'normal.bam')),
        mgd.TempOutputFile('tumour.unspiked.bam'),
        0.5)

    pyp.sch.commandline('simulate', (), ctx,
        os.path.join(tools_directory, 'bamextractsimreads'),
        '-b', mgd.InputFile(args['bam']),
        '-r', mgd.InputFile(args['ref']),
        '-s', mgd.InputFile(os.path.join(args['outdir'], 'simulated.fasta')),
        '-f', mgd.TempInputObj('simulation.params').extract(lambda a: a['coverage_fraction']),
        '-1', mgd.OutputFile(os.path.join(args['outdir'], 'simulated.1.fastq')),
        '-2', mgd.OutputFile(os.path.join(args['outdir'], 'simulated.2.fastq')))

    bwaalign_script = os.path.join(destruct_directory, 'bwaalign.py')

    pyp.sch.commandline('bwa_align', (), ctx,
        sys.executable,
        bwaalign_script,
        mgd.InputFile(args['ref']),
        mgd.InputFile(os.path.join(args['outdir'], 'simulated.1.fastq')),
        mgd.InputFile(os.path.join(args['outdir'], 'simulated.2.fastq')),
        mgd.TempOutputFile('simulated.unsorted.bam'),
        '--tmp', mgd.TempFile('bwa_tmp'))

    pyp.sch.transform('samtools_merge_sort_index', (), ctx,
        destruct_bam_test.samtools_merge_sort_index,
        None,
        mgd.OutputFile(os.path.join(args['outdir'], 'tumour.bam')),
        mgd.TempInputFile('tumour.unspiked.bam'),
        mgd.TempInputFile('simulated.unsorted.bam'))

    pyp.sch.transform('create_tool_wrappers', (), ctx,
        destruct_bam_test.create_tool_wrappers,
        mgd.TempOutputObj('tool_wrapper', 'bytool'),
        args['installdir'])

    pyp.sch.transform('run_tool', ('bytool',), ctx,
        destruct_bam_test.run_tool,
        None,
        mgd.TempInputObj('tool_wrapper', 'bytool'),
        mgd.TempFile('tool_tmp', 'bytool'),
        mgd.InputFile(os.path.join(args['outdir'], 'normal.bam')),
        mgd.InputFile(os.path.join(args['outdir'], 'tumour.bam')),
        mgd.OutputFile(os.path.join(args['outdir'], 'results_{bytool}.tsv'), 'bytool'))

    pyp.sch.transform('plot', ('bytool',), ctx, destruct_bam_test.create_roc_plot,
        None,
        mgd.TempInputObj('simulation.params'),
        mgd.TempInputObj('tool_wrapper', 'bytool'),
        mgd.InputFile(os.path.join(args['outdir'], 'simulated.tsv')),
        mgd.InputFile(os.path.join(args['outdir'], 'results_{bytool}.tsv'), 'bytool'),
        mgd.OutputFile(os.path.join(args['outdir'], 'annotated_{bytool}.tsv'), 'bytool'),
        mgd.OutputFile(os.path.join(args['outdir'], 'identified_{bytool}.tsv'), 'bytool'),
        mgd.OutputFile(os.path.join(args['outdir'], 'plots_{bytool}.pdf'), 'bytool'))

    pyp.run()


else:


    def partition_bam(original_filename, output_a_filename, output_b_filename, fraction_a):

        pypeliner.commandline.execute(os.path.join(tools_directory, 'bampartition'),
                                      '-i', original_filename,
                                      '-a', output_a_filename,
                                      '-b', output_b_filename,
                                      '-f', fraction_a)

        pypeliner.commandline.execute('samtools', 'index', output_a_filename, output_a_filename[:-4]+'.bai')
        pypeliner.commandline.execute('samtools', 'index', output_b_filename, output_b_filename[:-4]+'.bai')


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


    def create_tool_wrappers(install_directory):

        tool_wrappers = dict()

        for tool_name, ToolWrapper in wrappers.catalog.iteritems():

            tool_wrappers[tool_name] = ToolWrapper(os.path.join(install_directory, tool_name))

        return tool_wrappers


    def run_tool(tool_wrapper, temp_directory, normal_filename, tumour_filename, results_filename):

        tool_wrapper.run(temp_directory, {'normal':normal_filename, 'tumour':tumour_filename}, results_filename)


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

            results = results[results['normal_count'] == 0]

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

            num_missed = num_positive - np.sum(y_test)

            y_test = np.concatenate([y_test, np.array([True]*num_missed)])

            if y_test.sum() == len(y_test):
                y_test = np.concatenate([y_test, np.array([False])])

            fig = plt.figure(figsize=(16,16))

            for feature in features:
                values = results[feature].values
                values = np.concatenate([values, np.array([invalid_value] * (len(y_test) - len(values)))])
                fpr, tpr, thresholds = roc_curve(y_test, values)
                roc_auc = auc(fpr, tpr)
                plt.plot(fpr, tpr, label='{0} AUC = {1:.2f}'.format(feature, roc_auc))

            plt.plot([0, 1], [0, 1], 'k--')
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.0])
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title('Receiver operating characteristic example')
            plt.legend(loc="lower right")

            pdf.savefig(fig)

            for feature in features:

                fig = plt.figure(figsize=(16,16))

                true_features = results.loc[results['true_pos_id'].notnull(), feature].values.astype(float)
                false_features = results.loc[results['true_pos_id'].isnull(), feature].values.astype(float)

                def add_optional_noise(feature_values):
                    if len(np.unique(feature_values)) <= 1:
                        feature_values += np.random.randn(*feature_values.shape) * 1e-16
                    return feature_values

                true_features = add_optional_noise(true_features)
                false_features = add_optional_noise(false_features)

                seaborn.violinplot([true_features, false_features],
                                   names=['True', 'False'])

                plt.title('True vs false for ' + feature)

                pdf.savefig(fig)




