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
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.metrics import roc_curve, auc

import pypeliner

import demix

if __name__ == '__main__':

    import debalance_test
    import create_chromosome_simulation

    argparser = argparse.ArgumentParser()
    pypeliner.easypypeliner.add_arguments(argparser)
    argparser.add_argument('simconfig', help='Simulation configuration filename')
    argparser.add_argument('outdir', help='Output directory')
    argparser.add_argument('results', help='Test results plots (pdf)')

    cfg = pypeliner.easypypeliner.Config(vars(argparser.parse_args()))
    pyp = pypeliner.easypypeliner.EasyPypeliner([debalance_test, create_chromosome_simulation], cfg)

    if not cfg.results.lower().endswith('.pdf'):
        raise Exception('results file requires pdf extension')

    ctx = {'mem':4}

    try:
        os.makedirs(cfg.outdir)
    except OSError:
        pass

    pyp.sch.transform('read_params', (), ctx, debalance_test.read_simulation_params,
        pyp.sch.oobj('simulation.params'), 
        pyp.sch.input(cfg.simconfig))

    pyp.sch.transform('create_sim', (), ctx, create_chromosome_simulation.create, 
        None,
        cfg,
        pyp.sch.iobj('simulation.params'),
        pyp.sch.output(os.path.join(cfg.outdir, 'chromosome_info.tsv')),
        pyp.sch.output(os.path.join(cfg.outdir, 'breakpoint_info.tsv')),
        pyp.sch.output(os.path.join(cfg.outdir, 'simulated.1.fastq')),
        pyp.sch.output(os.path.join(cfg.outdir, 'simulated.2.fastq')),
        pyp.sch.tmpfile('simtemps'))

    bwaalign_script = os.path.join(cfg.destruct_directory, 'bwaalign.py')

    pyp.sch.commandline('bwa_align', (), ctx, 
        sys.executable,
        bwaalign_script,
        pyp.sch.input(os.path.join(cfg.outdir, 'simulated.1.fastq')),
        pyp.sch.input(os.path.join(cfg.outdir, 'simulated.2.fastq')),
        pyp.sch.output(os.path.join(cfg.outdir, 'simulated.bam')),  # This should be temp
        '--config', pyp.sch.input(cfg.config),
        '--tmp', pyp.sch.tmpfile('bwa_tmp'))

    pyp.sch.commandline('bam_sort', (), ctx,
        cfg.samtools_bin, 'sort', '-o',
        pyp.sch.input(os.path.join(cfg.outdir, 'simulated.bam')),
        pyp.sch.tmpfile('sorttemp'),
        '>', pyp.sch.output(os.path.join(cfg.outdir, 'simulated.sorted.bam')))

    pyp.sch.commandline('bam_index', (), ctx,
        cfg.samtools_bin, 'index',
        pyp.sch.input(os.path.join(cfg.outdir, 'simulated.sorted.bam')),
        pyp.sch.output(os.path.join(cfg.outdir, 'simulated.sorted.bam.bai')))

    pyp.sch.transform('write_bam_list', (), ctx, debalance_test.write_bam_list,
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

    pyp.sch.transform('create_changepoints', (), ctx, debalance_test.create_changepoints,
        None,
        pyp.sch.input(os.path.join(cfg.outdir, 'breakpoints.tsv')),
        pyp.sch.output(os.path.join(cfg.outdir, 'changepoints.tsv')))

    pyp.sch.commandline('bam_stats', (), ctx, cfg.bamstats_tool,
        '-b', pyp.sch.input(os.path.join(cfg.outdir, 'simulated.sorted.bam')),
        '--flen', '1000',
        '-s', pyp.sch.ofile('bamstats.file'))

    pyp.sch.transform('read_bam_stats', (), ctx, demix.read_stats,
        pyp.sch.oobj('bamstats'),
        pyp.sch.ifile('bamstats.file'))

    for chromosome in cfg.chromosomes.split():

        pyp.sch.commandline('read_concordant_{0}'.format(chromosome), (), ctx, cfg.bamconcordantreads_tool,
                            '--clipmax', '8', '--flen', '1000', '--chr', chromosome,
                            '-b', pyp.sch.input(os.path.join(cfg.outdir, 'simulated.sorted.bam')),
                            '-r', pyp.sch.ofile('reads.{0}'.format(chromosome)),
                            '-a', pyp.sch.ofile('alleles.{0}'.format(chromosome)))
        
        pyp.sch.transform('infer_haps_{0}'.format(chromosome), (), ctx, demix.infer_haps, None,
                          cfg, pyp.sch.temps_dir, pyp.sch.tmpfile('infer_haps'), chromosome, cfg.snp_positions,
                          pyp.sch.ifile('alleles.{0}'.format(chromosome)),
                          pyp.sch.ofile('hets.{0}'.format(chromosome)),
                          pyp.sch.ofile('haps.{0}'.format(chromosome)))

        pyp.sch.transform('create_readcounts_{0}'.format(chromosome), (), ctx, 
                          demix.create_counts, None, chromosome, 
                          pyp.sch.input(os.path.join(cfg.outdir, 'changepoints.tsv')),
                          pyp.sch.ifile('haps.{0}'.format(chromosome)),
                          pyp.sch.ifile('reads.{0}'.format(chromosome)),
                          pyp.sch.ifile('alleles.{0}'.format(chromosome)),
                          pyp.sch.ofile('interval.readcounts.{0}'.format(chromosome)),
                          pyp.sch.ofile('alleles.readcounts.{0}'.format(chromosome)),
                          pyp.sch.input(cfg.genome_fai))

    pyp.sch.transform('merge_interval_readcounts', (), ctx, demix.merge_files, None,
                      pyp.sch.ofile('interval.readcounts'),
                      *[pyp.sch.ifile('interval.readcounts.{0}'.format(chromosome)) for chromosome in cfg.chromosomes.split()])

    pyp.sch.commandline('samplegc', (), ctx, cfg.samplegc_tool, 
        '-b', pyp.sch.input(os.path.join(cfg.outdir, 'simulated.sorted.bam')),
        '-m', cfg.mappability_filename,
        '-g', cfg.genome_fasta,
        '-o', '4', '-n', '10000000',
        '-f', pyp.sch.iobj('bamstats').prop('fragment_length'),
        '>', pyp.sch.ofile('gcsamples'))

    pyp.sch.commandline('gcloess', (), ctx, cfg.rscript_bin, cfg.gc_loess_rscript,
        pyp.sch.ifile('gcsamples'),
        pyp.sch.ofile('gcloess'),
        pyp.sch.ofile('gcplots'))

    pyp.sch.commandline('gc_interval', (), ctx, cfg.estimategc_tool,
        '-m', cfg.mappability_filename,
        '-g', cfg.genome_fasta,
        '-c', pyp.sch.ifile('interval.readcounts'),
        '-i', '-o', '4', '-u',
        pyp.sch.iobj('bamstats').prop('fragment_mean'),
        '-s', pyp.sch.iobj('bamstats').prop('fragment_stddev'),
        '-a', cfg.mappability_length,
        '-l', pyp.sch.ifile('gcloess'),
        '>', pyp.sch.output(os.path.join(cfg.outdir, 'interval.readcounts.lengths')))



    # pyp.sch.transform('plot', (), ctx, debalance_test.create_roc_plot,
    #     None,
    #     pyp.sch.iobj('simulation.params'),
    #     pyp.sch.input(os.path.join(cfg.outdir, 'simulated.tsv')),
    #     pyp.sch.input(os.path.join(cfg.outdir, 'breakpoints.tsv')),
    #     pyp.sch.output(os.path.join(cfg.outdir, 'annotated.tsv')),
    #     pyp.sch.output(os.path.abspath(cfg.results)))

    pyp.run()

else:

    def read_simulation_params(config_filename):

        config = ConfigParser.ConfigParser()
        config.read(config_filename)
        return dict(config.items('main'))


    def write_bam_list(bam_list_filename, **kwargs):

        with open(bam_list_filename, 'w') as bam_list_file:
            for lib_id, bam_filename in kwargs.iteritems():
                bam_list_file.write(lib_id + '\t' + bam_filename + '\n')


    def create_changepoints(breakpoints_filename, changepoints_filename):

        breakpoints = pd.read_csv(breakpoints_filename, sep='\t', converters={'chromosome_1':str, 'chromosome_2':str})
        changepoints = pd.concat([breakpoints[['chromosome_1', 'break_1']].rename(columns=lambda a: a[:-2]), 
                                  breakpoints[['chromosome_2', 'break_2']].rename(columns=lambda a: a[:-2])], ignore_index=True)
        changepoints.to_csv(changepoints_filename, sep='\t', header=False, index=False)


    def create_roc_plot(sim_info, simulated_filename, predicted_filename, annotated_filename, plot_filename):

        breakpoints = pd.read_csv(simulated_filename, sep='\t', header=None,
                          converters={'chromosome1':str, 'chromosome2':str},
                          names=['break_id',
                                 'chromosome1', 'strand1', 'position1',
                                 'chromosome2', 'strand2', 'position2',
                                 'inserted', 'homology'])

        results = pd.read_csv(predicted_filename, sep='\t',
                              converters={'chromosome_1':str, 'chromosome_2':str})

        min_dist = 200

        def match(region, breakpoint):
            for side in (0, 1):
                if region[0][0] != breakpoint[side][0] or region[0][1] != breakpoint[side][1] or region[1][0] != breakpoint[1-side][0] or region[1][1] != breakpoint[1-side][1]:
                    continue
                if breakpoint[side][2] > region[0][3] + min_dist or breakpoint[side][2] < region[0][2] - min_dist or breakpoint[1-side][2] > region[1][3] + min_dist or breakpoint[1-side][2] < region[1][2] - min_dist:
                    continue
                return True
            return False

        def identify_true_positive(row):
            pred_region = [[row['chromosome_1'], row['strand_1'], row['start_1'], row['end_1']],
                           [row['chromosome_2'], row['strand_2'], row['start_2'], row['end_2']]]
            true_pos_id = None
            for idx, true_pos in breakpoints.iterrows():
                true_pos_region = [[true_pos['chromosome1'], true_pos['strand1'], true_pos['position1']],
                                   [true_pos['chromosome2'], true_pos['strand2'], true_pos['position2']]]
                if match(pred_region, true_pos_region):
                    true_pos_id = true_pos['break_id']
                    break
            return true_pos_id

        results['true_pos_id'] = results.apply(identify_true_positive, axis=1)

        results.to_csv(annotated_filename, sep='\t', index=False, na_rep='NA')

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

