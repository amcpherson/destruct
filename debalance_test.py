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
import breakpoint_graph
import eulerian


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

    pyp.sch.transform('bam_sort_index', (), ctx, debalance_test.bam_sort_index,
        None,
        cfg.samtools_bin,
        pyp.sch.input(os.path.join(cfg.outdir, 'simulated.bam')),
        pyp.sch.output(os.path.join(cfg.outdir, 'simulated.sorted.bam')),
        pyp.sch.output(os.path.join(cfg.outdir, 'simulated.sorted.bam.bai')),
        pyp.sch.tmpfile('sorttemp'))

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

    pyp.sch.transform('plot', (), ctx, debalance_test.infer_break_copies,
        None,
        pyp.sch.input(os.path.join(cfg.outdir, 'interval.readcounts.lengths')),
        pyp.sch.input(os.path.join(cfg.outdir, 'breakpoints.tsv')),
        pyp.sch.output(os.path.join(cfg.outdir, 'breakpoints_copies.tsv')))

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


    def bam_sort_index(samtools_bin, bam_filename, sorted_bam_filename, sorted_bam_index_filename, sort_temp):

        pypeliner.commandline.execute(samtools_bin, 'sort', '-o', bam_filename, sort_temp, '>', sorted_bam_filename)
        pypeliner.commandline.execute(samtools_bin, 'index', sorted_bam_filename, sorted_bam_index_filename)


    def create_changepoints(breakpoints_filename, changepoints_filename):

        breakpoints = pd.read_csv(breakpoints_filename, sep='\t', converters={'chromosome_1':str, 'chromosome_2':str})
        changepoints = pd.concat([breakpoints[['chromosome_1', 'break_1']].rename(columns=lambda a: a[:-2]), 
                                  breakpoints[['chromosome_2', 'break_2']].rename(columns=lambda a: a[:-2])], ignore_index=True)
        changepoints.to_csv(changepoints_filename, sep='\t', header=False, index=False)


    def infer_break_copies(intervals_filename, destruct_filename, results_filename):

        # Read interval data
        interval_data = pd.read_csv(intervals_filename, sep='\t', header=None, converters={'id':str, 'chromosome1':str, 'chromosome2':str},
                                    names=['id', 'chromosome1', 'position1', 'strand1', 'chromosome2', 'position2', 'strand2', 'readcount', 'length'])

        # Create reference data
        left_positions = interval_data[['chromosome1', 'position1']]
        right_positions = interval_data[['chromosome2', 'position2']]
        reference_data = pd.merge(left_positions, right_positions, left_on=['chromosome1', 'position1'], right_on=['chromosome2', 'position2'])
        reference_data['strand1'] = '+'
        reference_data['strand2'] = '-'
        reference_data['id'] = xrange(len(reference_data))

        # Read destruct results for variant data
        destruct_results = pd.read_csv(destruct_filename, sep='\t', converters={'chromosome_1':str, 'chromosome_2':str})
        variant_data = destruct_results[['cluster_id', 'chromosome_1', 'strand_1', 'break_1', 'chromosome_2', 'strand_2', 'break_2']]
        variant_data.columns = ['id', 'chromosome1', 'strand1', 'position1', 'chromosome2', 'strand2', 'position2']

        # Create telomeres as boundaries of intervals/variants
        def boundaries(df):
            return pd.concat([df[['chromosome1', 'strand1', 'position1']].rename(columns=lambda a: a[:-1]),
                              df[['chromosome2', 'strand2', 'position2']].rename(columns=lambda a: a[:-1])],
                             axis=0, ignore_index=True)

        interval_boundaries = boundaries(interval_data)

        variant_boundaries = boundaries(variant_data)
        variant_boundaries['is_breakpoint'] = True

        telomere_data = pd.merge(interval_boundaries, variant_boundaries, how='outer')
        telomere_data = telomere_data[pd.isnull(telomere_data['is_breakpoint'])]
        telomere_data = telomere_data.drop('is_breakpoint', axis=1)

        telomere_data = interval_boundaries

        telomere_data.set_index('chromosome', inplace=True)
        telomere_data['min_pos'] = telomere_data.groupby(level=0)['position'].min()
        telomere_data['max_pos'] = telomere_data.groupby(level=0)['position'].max()
        telomere_data.reset_index(inplace=True)
        telomere_data['is_normal'] = (telomere_data['position'] == telomere_data['min_pos']) | (telomere_data['position'] == telomere_data['max_pos'])
        telomere_data = telomere_data.drop(['min_pos', 'max_pos'], axis=1)
        telomere_data['id'] = xrange(len(telomere_data))

        bg = breakpoint_graph.BreakpointGraph(interval_data, reference_data, variant_data, telomere_data)

        ldas = np.linspace(0.0, 100000.0, num=101)

        coverage_paths = list()
        relevant = list()

        for lda in ldas:
            
            bg.clear_disabled_edges()

            solution = eulerian.solve_lasso(bg, lda, 1000.0 + 10.0 * lda)
            
            for edge_idx, coverage in enumerate(list(solution.edge_copies)):
                edge_type = bg.edges[edge_idx].edge_type
                if coverage <= 0.0 and (edge_type == 'variant' or edge_type == 'telomere'):
                    bg.set_disabled_edge(edge_idx)

            solution = eulerian.solve_lasso(bg, 0.0, 1000.0)

            if len(coverage_paths) == 0:
                num_edges = len(solution.edge_copies)
                coverage_paths = [list() for a in range(num_edges)]
                relevant = [True] * num_edges

            for edge_idx, coverage in enumerate(list(solution.edge_copies)):
                if coverage > 0.0:
                    if relevant[edge_idx]:
                        coverage_paths[edge_idx].append(coverage)
                else:
                    relevant[edge_idx] = False

        copies_results = list()

        hap_cov = interval_data['readcount'].sum() / interval_data['length'].sum()

        for edge_idx, coverage_path in enumerate(coverage_paths):

            if bg.edges[edge_idx].edge_type != 'variant':
                continue
            
            coverage_path = np.array(coverage_path)
            
            coverage_path = coverage_path[coverage_path > 0]
            
            relevance_score = float(len(coverage_path)) / float(len(ldas))
            
            if len(coverage_path) == 0:
                coverage_path = np.array([0])

            coverage_mean = np.mean(coverage_path)
            
            copies_mean = coverage_mean / hap_cov
            
            copies_results.append((bg.edges[edge_idx].id, coverage_mean, copies_mean, relevance_score))

        copies_results = pd.DataFrame(copies_results, columns=['cluster_id', 'coverage', 'copies', 'relevance'])

        copies_results = pd.merge(destruct_results, copies_results, left_on='cluster_id', right_on='cluster_id')

        copies_results.to_csv(results_filename, sep='\t', index=False)


