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

import pygenes
import pypeliner
import pypeliner.managed as mgd

import wrappers
import utils.download

destruct_directory = os.path.abspath(os.path.dirname(__file__))

default_config_filename = os.path.join(destruct_directory, 'defaultconfig.py')


if __name__ == '__main__':

    import destruct_sim_test
    import create_breakpoint_simulation

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('simconfig',
                           help='Simulation configuration filename')

    argparser.add_argument('installdir',
                           help='Tool installations directory')

    argparser.add_argument('outdir',
                           help='Output directory')

    argparser.add_argument('--config',
                           help='Configuration filename')

    argparser.add_argument('--chromosomes', nargs='*', type=str, default=['20'],
                           help='Reference chromosomes')

    argparser.add_argument('--include_nonchromosomal',  action='store_true',
                           help='Include non chromosomal reference sequences')

    args = vars(argparser.parse_args())

    config = {}

    if args['config'] is not None:
        execfile(args['config'], {}, config)

    config.update(args)

    pyp = pypeliner.app.Pypeline([destruct_sim_test, create_breakpoint_simulation], config)

    try:
        os.makedirs(args['outdir'])
    except OSError:
        pass

    ctx = {'mem':4}

    pyp.sch.transform('read_params', (), ctx, destruct_sim_test.read_simulation_params,
        mgd.TempOutputObj('simulation.params'),
        mgd.InputFile(args['simconfig']))

    pyp.sch.transform('create_genome', (), ctx,
        destruct_sim_test.create_genome,
        None,
        args['chromosomes'],
        args['include_nonchromosomal'],
        mgd.OutputFile(os.path.join(args['outdir'], 'genome.fasta')))

    pyp.sch.transform('create_sim', (), ctx, create_breakpoint_simulation.create,
        None,
        mgd.TempInputObj('simulation.params'),
        mgd.InputFile(os.path.join(args['outdir'], 'genome.fasta')),
        mgd.OutputFile(os.path.join(args['outdir'], 'simulated.fasta')),
        mgd.OutputFile(os.path.join(args['outdir'], 'simulated.tsv')),
        mgd.TempOutputFile('concordant.1.fastq'),
        mgd.TempOutputFile('concordant.2.fastq'),
        mgd.TempOutputFile('discordant.1.fastq'),
        mgd.TempOutputFile('discordant.2.fastq'))

    pyp.sch.commandline('cat1', (), ctx, 'cat', mgd.TempInputFile('concordant.1.fastq'), mgd.TempInputFile('discordant.1.fastq'), '>', mgd.OutputFile(os.path.join(args['outdir'], 'simulated.1.fastq')))
    pyp.sch.commandline('cat2', (), ctx, 'cat', mgd.TempInputFile('concordant.2.fastq'), mgd.TempInputFile('discordant.2.fastq'), '>', mgd.OutputFile(os.path.join(args['outdir'], 'simulated.2.fastq')))

    bwaalign_script = os.path.join(destruct_directory, 'bwaalign.py')

    pyp.sch.commandline('bwa_align', (), ctx, 
        sys.executable,
        bwaalign_script,
        mgd.InputFile(os.path.join(args['outdir'], 'genome.fasta')),
        mgd.InputFile(os.path.join(args['outdir'], 'simulated.1.fastq')),
        mgd.InputFile(os.path.join(args['outdir'], 'simulated.2.fastq')),
        mgd.TempOutputFile('simulated.unsorted.bam'),
        '--tmp', mgd.TempFile('bwa_tmp'))

    pyp.sch.transform('samtools_sort_index', (), ctx,
        destruct_sim_test.samtools_sort_index,
        None,
        mgd.TempInputFile('simulated.unsorted.bam'),
        mgd.OutputFile(os.path.join(args['outdir'], 'simulated.bam')))

    pyp.sch.transform('set_tools', (), ctx,
        destruct_sim_test.set_tools,
        mgd.OutputChunks('bytool'))

    pyp.sch.transform('run_tool', ('bytool',), ctx,
        destruct_sim_test.run_tool,
        None,
        mgd.InputInstance('bytool'),
        mgd.Template(os.path.join(args['installdir'], '{bytool}'), 'bytool'),
        mgd.TempFile('tool_tmp', 'bytool'),
        mgd.InputFile(os.path.join(args['outdir'], 'simulated.bam')),
        mgd.OutputFile(os.path.join(args['outdir'], 'results_{bytool}.tsv'), 'bytool'))

    pyp.sch.transform('plot', ('bytool',), ctx, destruct_sim_test.create_roc_plot,
        None,
        mgd.TempInputObj('simulation.params'),
        mgd.InputFile(os.path.join(args['outdir'], 'simulated.tsv')),
        mgd.InputFile(os.path.join(args['outdir'], 'results_{bytool}.tsv'), 'bytool'),
        mgd.OutputFile(os.path.join(args['outdir'], 'annotated_{bytool}.tsv'), 'bytool'),
        mgd.OutputFile(os.path.join(args['outdir'], 'identified_{bytool}.tsv'), 'bytool'),
        mgd.OutputFile(os.path.join(args['outdir'], 'plots_{bytool}.pdf'), 'bytool'))

    pyp.run()


else:


    def create_genome(chromosomes, include_nonchromosomal, genome_fasta):

        utils.download.download_genome_fasta(genome_fasta,
                                             chromosomes,
                                             include_nonchromosomal)

        pypeliner.commandline.execute('bwa', 'index', genome_fasta)
        pypeliner.commandline.execute('samtools', 'faidx', genome_fasta)

        for extension in ('.fai', '.amb', '.ann', '.bwt', '.pac', '.sa'):
            os.rename(genome_fasta+extension, genome_fasta[:-4]+extension)


    def samtools_sort_index(input_filename, output_filename):

        pypeliner.commandline.execute('samtools', 'sort', input_filename, output_filename)

        os.rename(output_filename + '.bam', output_filename)

        assert output_filename.endswith('.tmp')
        index_filename = output_filename[:-4] + '.bai'

        pypeliner.commandline.execute('samtools', 'index', output_filename, index_filename)


    def set_tools():
        return wrappers.catalog.keys()


    def run_tool(tool, install_directory, temp_directory, bam_filename, results_filename):

        try:
            ToolWrapper = wrappers.catalog[tool]
        except KeyError:
            raise Exception('No wrapper for tool ' + tool)

        tool_wrapper = ToolWrapper(install_directory)

        bam_filenames = {'simulated':bam_filename}

        tool_wrapper.run(bam_filenames, results_filename, temp_directory)


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


    def create_roc_plot(sim_info, simulated_filename, predicted_filename, annotated_filename, identified_filename, plot_filename):

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

        features = ['simulated_count', 'num_split']
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

        fig.savefig(plot_filename, format='pdf')

