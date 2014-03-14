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

if __name__ == '__main__':

    import destruct_test
    import create_breakpoint_simulation

    argparser = argparse.ArgumentParser()
    pypeliner.easypypeliner.add_arguments(argparser)
    argparser.add_argument('simconfig', help='Simulation configuration filename')
    argparser.add_argument('outdir', help='Output directory')
    argparser.add_argument('results', help='Test results plots (pdf)')

    cfg = pypeliner.easypypeliner.Config(vars(argparser.parse_args()))
    pyp = pypeliner.easypypeliner.EasyPypeliner([destruct_test, create_breakpoint_simulation], cfg)

    if not cfg.results.lower().endswith('.pdf'):
        raise Exception('results file requires pdf extension')

    ctx = {'mem':4}

    try:
        os.makedirs(cfg.outdir)
    except OSError:
        pass

    pyp.sch.transform('read_params', (), ctx, destruct_test.read_simulation_params,
        pyp.sch.oobj('simulation.params'), 
        pyp.sch.input(cfg.simconfig))

    pyp.sch.transform('create_sim', (), ctx, create_breakpoint_simulation.create, 
        None,
        cfg,
        pyp.sch.iobj('simulation.params'),
        pyp.sch.output(os.path.join(cfg.outdir, 'simulated.fasta')),
        pyp.sch.output(os.path.join(cfg.outdir, 'simulated.tsv')),
        pyp.sch.ofile('concordant.1.fastq'),
        pyp.sch.ofile('concordant.2.fastq'),
        pyp.sch.ofile('discordant.1.fastq'),
        pyp.sch.ofile('discordant.2.fastq'))

    pyp.sch.commandline('cat1', (), ctx, 'cat', pyp.sch.ifile('concordant.1.fastq'), pyp.sch.ifile('discordant.1.fastq'), '>', pyp.sch.output(os.path.join(cfg.outdir, 'simulated.1.fastq')))
    pyp.sch.commandline('cat2', (), ctx, 'cat', pyp.sch.ifile('concordant.2.fastq'), pyp.sch.ifile('discordant.2.fastq'), '>', pyp.sch.output(os.path.join(cfg.outdir, 'simulated.2.fastq')))

    bwaalign_script = os.path.join(cfg.destruct_directory, 'bwaalign.py')

    pyp.sch.commandline('bwa_align', (), ctx, 
        sys.executable,
        bwaalign_script,
        pyp.sch.input(os.path.join(cfg.outdir, 'simulated.1.fastq')),
        pyp.sch.input(os.path.join(cfg.outdir, 'simulated.2.fastq')),
        pyp.sch.output(os.path.join(cfg.outdir, 'simulated.bam')),
        '--config', pyp.sch.input(cfg.config),
        '--tmp', pyp.sch.tmpfile('bwa_tmp'))

    pyp.sch.transform('write_bam_list', (), ctx, destruct_test.write_bam_list,
        None,
        pyp.sch.output(os.path.join(cfg.outdir, 'bam_list.tsv')),
        **{'simulated':pyp.sch.input(os.path.join(cfg.outdir, 'simulated.bam'))})

    destruct_script = os.path.join(cfg.destruct_directory, 'destruct.py')

    pyp.sch.commandline('destruct', (), ctx, 
        sys.executable,
        destruct_script,
        pyp.sch.input(os.path.join(cfg.outdir, 'bam_list.tsv')),
        pyp.sch.output(os.path.join(cfg.outdir, 'breakpoints.tsv')),
        pyp.sch.output(os.path.join(cfg.outdir, 'breakreads.tsv')),
        '--config', pyp.sch.input(cfg.config),
        '--tmp', pyp.sch.tmpfile('destruct_tmp'),
        '--nocleanup')

    pyp.sch.transform('plot', (), ctx, destruct_test.create_roc_plot,
        None,
        pyp.sch.iobj('simulation.params'),
        pyp.sch.input(os.path.join(cfg.outdir, 'simulated.tsv')),
        pyp.sch.input(os.path.join(cfg.outdir, 'breakpoints.tsv')),
        pyp.sch.output(os.path.join(cfg.outdir, 'annotated.tsv')),
        pyp.sch.output(os.path.join(cfg.outdir, 'identified.tsv')),
        pyp.sch.output(os.path.abspath(cfg.results)))

    pyp.run()

else:

    def read_simulation_params(config_filename):

        config = ConfigParser.ConfigParser()
        config.read(config_filename)
        section = 'main'

        sim_info = dict()
        sim_info['read_length'] = config.get(section, 'read_length')
        sim_info['fragment_mean'] = config.get(section, 'fragment_mean')
        sim_info['fragment_stddev'] = config.get(section, 'fragment_stddev')
        sim_info['coverage'] = config.get(section, 'coverage')
        sim_info['num_inserted'] = config.get(section, 'num_inserted')
        sim_info['homology'] = config.get(section, 'homology')
        sim_info['adjacent_length'] = config.get(section, 'adjacent_length')
        sim_info['breakpoints_seed'] = config.get(section, 'breakpoints_seed')
        sim_info['dwgsim_seed'] = config.get(section, 'dwgsim_seed')
        sim_info['num_breakpoints'] = config.get(section, 'num_breakpoints')
        sim_info['num_concordant'] = config.get(section, 'num_concordant')

        return sim_info


    def write_bam_list(bam_list_filename, **kwargs):

        with open(bam_list_filename, 'w') as bam_list_file:
            for lib_id, bam_filename in kwargs.iteritems():
                bam_list_file.write(lib_id + '\t' + bam_filename + '\n')


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
                idx = bisect.bisect_left(chrom_strand_positions, row['break_'+side] - extend)
                side_matched_ids = list()
                while idx < len(chrom_strand_positions):
                    pos = chrom_strand_positions[idx]
                    dist = abs(pos - row['break_'+side])
                    if pos >= row['break_'+side] - extend and pos <= row['break_'+side] + extend:
                        for break_id in self.break_ids[(row['chromosome_'+side], row['strand_'+side], pos)]:
                            side_matched_ids.append((break_id, dist))
                    if pos > row['break_'+side] + extend:
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
                              converters={'cluster_id':str, 'chromosome_1':str, 'chromosome_2':str})

        min_dist = 200

        def identify_true_positive(row):
            return breakpoints_db.query(row, min_dist)

        results['true_pos_id'] = results.apply(identify_true_positive, axis=1)

        results.to_csv(annotated_filename, sep='\t', index=False, na_rep='NA')

        identified_columns = ['align_prob', 'valid_prob', 'chimeric_prob', 'simulated_count', 'num_split']

        identified = breakpoints[['break_id']]
        identified = identified.merge(results[['cluster_id', 'true_pos_id'] + identified_columns], left_on='break_id', right_on='true_pos_id', how='outer')
        identified = identified.drop('true_pos_id', axis=1)

        identified.to_csv(identified_filename, sep='\t', index=False, na_rep='NA')

        num_positive = int(sim_info['num_breakpoints'])

        y_test = results['true_pos_id'].notnull()

        num_missed = num_positive - np.sum(y_test)

        y_test = np.concatenate([y_test, np.array([False]*num_missed)])

        fig = plt.figure(figsize=(16,16))

        for feature in ('align_prob', 'valid_prob', 'chimeric_prob', 'simulated_count', 'num_split'):
            probs = np.concatenate([results[feature], np.array([0.0]*num_missed)])
            fpr, tpr, thresholds = roc_curve(y_test, probs)
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

