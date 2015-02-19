import csv
import sys
import logging
import os
import ConfigParser
import itertools
import argparse
import string
import gzip
import shutil
import tarfile
import io
from collections import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import math
import scipy.stats
import sklearn
import sklearn.mixture
import statsmodels.api as sm

import pypeliner
import pypeliner.managed as mgd

import subclonal_sampling


destruct_directory = os.path.abspath(os.path.dirname(__file__))

data_directory = os.path.join(destruct_directory, 'data')
tools_directory = os.path.join(destruct_directory, 'tools')
default_config_filename = os.path.join(destruct_directory, 'defaultconfig.py')


__version__ = '0.0.1'

if __name__ == '__main__':

    import demix
    
    argparser = argparse.ArgumentParser()

    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('--version', action='version', version=__version__)

    argparser.add_argument('ref_data_dir',
                           help='Reference dataset directory')

    argparser.add_argument('libraries',
                           help='library info')

    argparser.add_argument('normal',
                           help='normal library id')

    argparser.add_argument('changepoints',
                           help='changepoints file')

    argparser.add_argument('stats',
                           help='statistics file')

    argparser.add_argument('preds',
                           help='predictions file')

    argparser.add_argument('plots',
                           help='plots file tarball')

    argparser.add_argument('plots_prefix',
                           help='plots filename prefix')

    argparser.add_argument('--config', required=False,
                           help='Configuration Filename')

    args = vars(argparser.parse_args())

    config = {'ref_data_directory':args['ref_data_dir'],
              'package_data_directory':data_directory}
    execfile(default_config_filename, {}, config)

    if args['config'] is not None:
        execfile(args['config'], {}, config)

    config.update(args)

    pyp = pypeliner.app.Pypeline([demix], config)

    ctx_general = {'mem':16, 'ncpus':1}

    tumour_lib_infos = demix.read_libs(config['libraries'])

    try:
        normal_lib_info = tumour_lib_infos[config['normal']]
        del tumour_lib_infos[config['normal']]
    except KeyError:
        raise Exception('normal library id "{0}" should match one entry in library list "{1}"'.format(config['normal'], config['libraries']))

    for lib_info in tumour_lib_infos.values() + [normal_lib_info]:

        pyp.sch.commandline('bam_stats_'+lib_info.id, (), ctx_general,
            os.path.join(tools_directory, 'bamstats'),
            '-b', mgd.InputFile(lib_info.bam_filename),
            '--flen', '1000',
            '-s', mgd.TempOutputFile('bamstats.file.'+lib_info.id))

        pyp.sch.transform('read_bam_stats_'+lib_info.id, (), ctx_general,
            demix.read_stats,
            mgd.TempOutputObj('bamstats.'+lib_info.id),
            mgd.TempInputFile('bamstats.file.'+lib_info.id))

    for chromosome in config['chromosomes']:

        pyp.sch.commandline('read_concordant_{0}_{1}'.format(chromosome, normal_lib_info.id), (), ctx_general,
            os.path.join(tools_directory, 'bamconcordantreads'),
            '--clipmax', '8',
            '--flen', '1000',
            '--chr', chromosome,
            '-b', mgd.InputFile(normal_lib_info.bam_filename),
            '-s', config['snp_positions'],
            '-r', mgd.TempOutputFile('reads.{0}.{1}'.format(chromosome, normal_lib_info.id)),
            '-a', mgd.TempOutputFile('alleles.{0}.{1}'.format(chromosome, normal_lib_info.id)))

        pyp.sch.transform('infer_haps_{0}'.format(chromosome), (), ctx_general,
            demix.infer_haps,
            None,
            config,
            pyp.sch.temps_dir,
            normal_lib_info.id,
            chromosome,
            config['snp_positions'],
            mgd.TempInputFile('alleles.{0}.{1}'.format(chromosome, normal_lib_info.id)),
            mgd.TempOutputFile('hets.{0}'.format(chromosome)),
            mgd.TempOutputFile('haps.{0}'.format(chromosome)))

    for lib_info in tumour_lib_infos.values():

        for chromosome in config['chromosomes']:

            pyp.sch.commandline('read_concordant_{0}_{1}'.format(chromosome, lib_info.id), (), ctx_general,
                os.path.join(tools_directory, 'bamconcordantreads'),
                '--clipmax', '8',
                '--flen', '1000',
                '--chr', chromosome,
                '-b', mgd.InputFile(lib_info.bam_filename),
                '-s', mgd.TempInputFile('hets.{0}'.format(chromosome)),
                '-r', mgd.TempOutputFile('reads.{0}.{1}'.format(chromosome, lib_info.id)),
                '-a', mgd.TempOutputFile('alleles.{0}.{1}'.format(chromosome, lib_info.id)))

    for lib_info in tumour_lib_infos.values():

        for chromosome in config['chromosomes']:

            pyp.sch.transform('create_readcounts_{0}_{1}'.format(chromosome, lib_info.id), (), ctx_general,
                demix.create_counts,
                None,
                chromosome, 
                mgd.InputFile(config['changepoints']),
                mgd.TempInputFile('haps.{0}'.format(chromosome)),
                mgd.TempInputFile('reads.{0}.{1}'.format(chromosome, lib_info.id)),
                mgd.TempInputFile('alleles.{0}.{1}'.format(chromosome, lib_info.id)),
                mgd.TempOutputFile('interval.readcounts.{0}.{1}'.format(chromosome, lib_info.id)),
                mgd.TempOutputFile('alleles.readcounts.{0}.{1}'.format(chromosome, lib_info.id)),
                mgd.InputFile(config['genome_fai']))

    for chromosome in config['chromosomes']:

        pyp.sch.transform('phase_intervals_{0}'.format(chromosome), (), ctx_general,
            demix.phase_intervals,
            None,
            *([mgd.TempInputFile('alleles.readcounts.{0}.{1}'.format(chromosome, lib_info.id)) for lib_info in sorted(tumour_lib_infos.values())] + 
              [mgd.TempOutputFile('alleles.readcounts.phased.{0}.{1}'.format(chromosome, lib_info.id)) for lib_info in sorted(tumour_lib_infos.values())]))

    for lib_info in tumour_lib_infos.values():

        pyp.sch.transform('merge_interval_readcounts_{0}'.format(lib_info.id), (), ctx_general,
            demix.merge_files,
            None,
            mgd.TempOutputFile('interval.readcounts.{0}'.format(lib_info.id)),
            *[mgd.TempInputFile('interval.readcounts.{0}.{1}'.format(chromosome, lib_info.id)) for chromosome in config['chromosomes']])

        pyp.sch.transform('merge_allele_readcounts_{0}'.format(lib_info.id), (), ctx_general,
            demix.merge_files,
            None,
            mgd.TempOutputFile('alleles.readcounts.phased.{0}'.format(lib_info.id)),
            *[mgd.TempInputFile('alleles.readcounts.phased.{0}.{1}'.format(chromosome, lib_info.id)) for chromosome in config['chromosomes']])

        pyp.sch.commandline('samplegc_'+lib_info.id, (), ctx_general,
            os.path.join(tools_directory, 'samplegc'),
            '-b', mgd.InputFile(lib_info.bam_filename),
            '-m', config['mappability_filename'],
            '-g', config['genome_fasta'],
            '-o', '4',
            '-n', '10000000',
            '-f', mgd.TempInputObj('bamstats.'+lib_info.id).prop('fragment_length'),
            '>', mgd.TempOutputFile('gcsamples.'+lib_info.id))

        pyp.sch.transform('gcloess_'+lib_info.id, (), ctx_general,
            demix.gc_lowess,
            None,
            mgd.TempInputFile('gcsamples.'+lib_info.id),
            mgd.TempOutputFile('gcloess.'+lib_info.id),
            mgd.TempOutputFile('gcplots.'+lib_info.id))

        pyp.sch.commandline('gc_interval_'+lib_info.id, (), ctx_general,
            os.path.join(tools_directory, 'estimategc'),
            '-m', config['mappability_filename'],
            '-g', config['genome_fasta'],
            '-c', mgd.TempInputFile('interval.readcounts.'+lib_info.id),
            '-i',
            '-o', '4',
            '-u', mgd.TempInputObj('bamstats.'+lib_info.id).prop('fragment_mean'),
            '-s', mgd.TempInputObj('bamstats.'+lib_info.id).prop('fragment_stddev'),
            '-a', config['mappability_length'],
            '-l', mgd.TempInputFile('gcloess.'+lib_info.id),
            '>', mgd.TempOutputFile('interval.readcounts.lengths.'+lib_info.id))

        pyp.sch.transform('solve_and_plot_'+lib_info.id, (), ctx_general,
            demix.solve_and_plot,
            None,
            lib_info.id,
            mgd.TempInputFile('interval.readcounts.lengths.'+lib_info.id),
            mgd.TempInputFile('alleles.readcounts.phased.{0}'.format(lib_info.id)),
            mgd.TempOutputFile('stats.{0}'.format(lib_info.id)),
            mgd.TempOutputFile('preds.{0}'.format(lib_info.id)),
            mgd.TempOutputFile('plots.{0}'.format(lib_info.id)),
            config['plots_prefix'])

    pyp.sch.transform('merge_stats', (), ctx_general, demix.merge_tables, None, mgd.OutputFile(config['stats']),
        *[mgd.TempInputFile('stats.{0}'.format(lib_info.id)) for lib_info in tumour_lib_infos.values()])
    pyp.sch.transform('merge_preds', (), ctx_general, demix.merge_tables, None, mgd.OutputFile(config['preds']),
        *[mgd.TempInputFile('preds.{0}'.format(lib_info.id)) for lib_info in tumour_lib_infos.values()])
    pyp.sch.transform('merge_plots', (), ctx_general, demix.merge_tars, None, mgd.OutputFile(config['plots']),
        *[mgd.TempInputFile('plots.{0}'.format(lib_info.id)) for lib_info in tumour_lib_infos.values()])

    pyp.run()


LibInfo = namedtuple('LibInfo', ['id', 'bam_filename'])

def read_libs(libs_filename):
    lib_infos = dict()
    with open(libs_filename, 'r') as libs_file:
        for row in csv.reader(libs_file, delimiter='\t'):
            id = row[0]
            bam_filename = row[1]
            lib_infos[id] = LibInfo(id, bam_filename)
    return lib_infos


class ConcordantReadStats(object):
    def __init__(self, stats):
        self.stats = stats
    @property
    def fragment_mean(self):
        return float(self.stats['fragment_mean'])
    @property
    def fragment_stddev(self):
        return float(self.stats['fragment_stddev'])
    @property
    def fragment_length(self):
        return int(float(self.stats['fragment_mean']))


def read_stats(stats_filename):
    with open(stats_filename, 'r') as stats_file:
        header = stats_file.readline().rstrip().split('\t')
        values = stats_file.readline().rstrip().split('\t')
        return ConcordantReadStats(dict(zip(header,values)))


def is_contained(a, b):
    """ Check if region b is fully contained within region a """
    return b[0] >= a[0] and b[1] <= a[1]


def contained_counts(X, Y):
    """ Find counts of regions in Y contained in regions in X
    X and Y are assumed sorted by start
    X is a set of non-overlapping regions
    """
    C = np.zeros(X.shape[0])
    y_idx = 0
    for x_idx, x in enumerate(X):
        while y_idx < Y.shape[0] and Y[y_idx][0] < x[0]:
            y_idx += 1
        while y_idx < Y.shape[0] and Y[y_idx][0] <= x[1]:
            if is_contained(x, Y[y_idx]):
                C[x_idx] += 1
            y_idx += 1
    return C


def overlapping_counts(X, Y):
    """ Find counts of regions in Y overlapping positions in X
    X and Y are assume sorted, Y by starting position, X by position
    """
    C = np.zeros(X.shape[0])
    x_idx = 0
    for y in Y:
        while x_idx < X.shape[0] and X[x_idx] <= y[0]:
            x_idx += 1
        x_idx_1 = x_idx
        while x_idx_1 < X.shape[0] and X[x_idx_1] < y[1]:
            C[x_idx_1] += 1
            x_idx_1 += 1
    return C


def find_contained(X, Y):
    """ Find mapping of positions in Y contained within regions in X
    X and Y are assume sorted, X by starting position, Y by position
    X is a set of non-overlapping regions
    """
    M = [None]*Y.shape[0]
    y_idx = 0
    for x_idx, x in enumerate(X):
        while y_idx < Y.shape[0] and Y[y_idx] <= x[1]:
            if Y[y_idx] >= x[0]:
                M[y_idx] = x_idx
            y_idx += 1
    return M


def read_reads_data(reads_filename, num_rows=-1):
    dt = np.dtype([('start', np.uint32), ('length', np.uint16)])
    with gzip.open(reads_filename, 'rb') as reads_file:
        read_id_start = 0
        while True:
            raw_data = reads_file.read(num_rows * dt.itemsize)
            if raw_data == '':
                yield pd.DataFrame(columns=['id', 'start', 'end'])
                break
            data = np.fromstring(raw_data, dtype=dt)
            df = pd.DataFrame(data)
            df['id'] = xrange(read_id_start, read_id_start+len(df))
            df['end'] = df['start'] + df['length'] - 1
            df = df.drop('length', axis=1)
            yield df


def read_alleles_data(alleles_filename, num_rows=-1):
    dt = np.dtype([('read', np.uint32), ('pos', np.uint32), ('is_alt', np.uint8)])
    with gzip.open(alleles_filename, 'rb') as alleles_file:
        while True:
            raw_data = alleles_file.read(num_rows * dt.itemsize)
            if raw_data == '':
                yield pd.DataFrame(columns=['read', 'pos', 'is_alt'])
                break
            data = np.fromstring(raw_data, dtype=dt)
            df = pd.DataFrame(data)
            yield df


def create_counts(chromosome, changepoints_filename, haps_filename, reads_filename, 
                  alleles_filename, interval_filename, allele_counts_filename,
                  genome_fai_filename):
    
    # Read changepoint data
    changepoints = pd.read_csv(changepoints_filename, sep='\t', header=None,
                               converters={'chromosome':str}, 
                               names=['chromosome', 'position'])
    haps = pd.read_csv(haps_filename, sep='\t')
    reads = next(read_reads_data(reads_filename))
    reads.sort('start', inplace=True)
    
    # Create a list of regions between changepoints
    changepoints = changepoints[changepoints['chromosome'] == chromosome]
    changepoints = changepoints.append({'chromosome':chromosome, 'position':1}, ignore_index=True)
    changepoints = changepoints.append({'chromosome':chromosome, 'position':read_chromosome_lengths(genome_fai_filename)[chromosome]}, ignore_index=True)
    changepoints.drop_duplicates(inplace=True)
    changepoints.sort('position', inplace=True)
    regions = pd.DataFrame(data=np.array([changepoints['position'].values[:-1],
                                          changepoints['position'].values[1:]]).T,
                           columns=['start', 'end'])
    regions.drop_duplicates(inplace=True)
    regions.sort('start', inplace=True)
    
    # Create an index that matches the sort order
    regions.index = xrange(len(regions))

     # Count interval reads
    interval_counts = contained_counts(regions[['start', 'end']].values, reads[['start', 'end']].values)

    del reads

    # Create interval data
    interval_data = pd.DataFrame({'start':regions['start'].values, 'end':regions['end'].values, 'counts':interval_counts})

    # Write interval data to a file
    interval_data['id'] = chromosome + '_'
    interval_data['id'] += interval_data.index.values.astype(str)
    interval_data['counts'] = interval_data['counts'].astype(int)
    interval_data['chromosome1'] = chromosome
    interval_data['strand1'] = '-'
    interval_data['chromosome2'] = chromosome
    interval_data['strand2'] = '+'
    interval_data = interval_data[['id', 'chromosome1', 'start', 'strand1', 'chromosome2', 'end', 'strand2', 'counts']]
    interval_data.to_csv(interval_filename, sep='\t', index=False, header=False)

    # Merge haplotype information into read alleles table
    alleles = list()
    for alleles_chunk in read_alleles_data(alleles_filename, num_rows=10000):
        alleles_chunk = alleles_chunk.merge(haps, left_on=['pos', 'is_alt'], right_on=['pos', 'allele'], how='inner')
        alleles.append(alleles_chunk)
    alleles = pd.concat(alleles, ignore_index=True)

    # Arbitrarily assign a haplotype/allele label to each read
    alleles.drop_duplicates('read', inplace=True)

    # Create a mapping between regions and snp positions
    snp_region = pd.DataFrame({'pos':haps['pos'].unique()})
    snp_region['region_idx'] = find_contained(regions[['start', 'end']].values, snp_region['pos'].values)
    snp_region = snp_region.dropna()
    snp_region['region_idx'] = snp_region['region_idx'].astype(int)

    # Add annotation of which region each snp is contained within
    alleles = alleles.merge(snp_region, left_on='pos', right_on='pos')

    # Count reads for each allele
    alleles.set_index(['region_idx', 'hap_label', 'allele_id'], inplace=True)
    allele_counts = alleles.groupby(level=[0, 1, 2]).size().reset_index().rename(columns={0:'count'})

    # Create region id as chromosome _ index
    allele_counts['region_id'] = chromosome + '_'
    allele_counts['region_id'] += allele_counts['region_idx'].astype(str)

    # Write out allele counts
    allele_counts.to_csv(allele_counts_filename, sep='\t', cols=['region_id', 'hap_label', 'allele_id', 'count'], index=False, header=False)


def infer_haps(config, temps_directory, library, chromosome, snps_filename, normal_alleles_filename, hets_filename, haps_filename):
    
    def write_null():
        with open(hets_filename, 'w') as hets_file:
            pass
        with open(haps_filename, 'w') as haps_file:
            haps_file.write('pos\tallele\tchangepoint_confidence\thap_label\tallele_id\tallele_label\n')

    accepted_chromosomes = [str(a) for a in range(1, 23)] + ['X']
    if str(chromosome) not in accepted_chromosomes:
        write_null()
        return
    
    # Temporary directory for impute2 files
    haps_temp_directory = os.path.join(os.path.join(temps_directory, library), chromosome)
    try:
        os.makedirs(haps_temp_directory)
    except OSError:
        pass

    # Impute 2 files for thousand genomes data by chromosome
    phased_chromosome = chromosome
    if chromosome == 'X':
        phased_chromosome = config['phased_chromosome_x']
    genetic_map_filename = config['genetic_map_template'].format(phased_chromosome)
    hap_filename = config['haplotypes_template'].format(phased_chromosome)
    legend_filename = config['legend_template'].format(phased_chromosome)

    # Call snps based on reference and alternate read counts from normal
    snp_counts_df = list()
    for alleles_chunk in read_alleles_data(normal_alleles_filename, num_rows=10000):
        snp_counts_chunk = alleles_chunk.groupby(['pos', 'is_alt']).size().unstack().fillna(0)
        snp_counts_chunk = snp_counts_chunk.rename(columns=lambda a: {0:'ref_count', 1:'alt_count'}[a])
        snp_counts_chunk = snp_counts_chunk.astype(float)
        snp_counts_df.append(snp_counts_chunk)
    snp_counts_df = pd.concat(snp_counts_df)
    snp_counts_df = snp_counts_df.groupby(level=0).sum()
    snp_counts_df.sort_index(inplace=True)

    if len(snp_counts_df) == 0:
        write_null()
        return

    snp_counts_df['total_count'] = snp_counts_df['ref_count'] + snp_counts_df['alt_count']

    snp_counts_df['likelihood_AA'] = scipy.stats.binom.pmf(snp_counts_df['alt_count'], snp_counts_df['total_count'], float(config['sequencing_base_call_error']))
    snp_counts_df['likelihood_AB'] = scipy.stats.binom.pmf(snp_counts_df['alt_count'], snp_counts_df['total_count'], 0.5)
    snp_counts_df['likelihood_BB'] = scipy.stats.binom.pmf(snp_counts_df['ref_count'], snp_counts_df['total_count'], float(config['sequencing_base_call_error']))
    snp_counts_df['evidence'] = snp_counts_df['likelihood_AA'] + snp_counts_df['likelihood_AB'] + snp_counts_df['likelihood_BB']

    snp_counts_df['posterior_AA'] = snp_counts_df['likelihood_AA'] / snp_counts_df['evidence']
    snp_counts_df['posterior_AB'] = snp_counts_df['likelihood_AB'] / snp_counts_df['evidence']
    snp_counts_df['posterior_BB'] = snp_counts_df['likelihood_BB'] / snp_counts_df['evidence']

    snp_counts_df['AA'] = (snp_counts_df['posterior_AA'] >= float(config['het_snp_call_threshold'])) * 1
    snp_counts_df['AB'] = (snp_counts_df['posterior_AB'] >= float(config['het_snp_call_threshold'])) * 1
    snp_counts_df['BB'] = (snp_counts_df['posterior_BB'] >= float(config['het_snp_call_threshold'])) * 1

    snp_counts_df = snp_counts_df[(snp_counts_df['AA'] == 1) | (snp_counts_df['AB'] == 1) | (snp_counts_df['BB'] == 1)]

    snps_df_iter = pd.read_csv(snps_filename, sep='\t', names=['chr', 'pos', 'ref', 'alt'], converters={'chr':str}, iterator=True, chunksize=10000)
    snps_df = pd.concat([chunk[chunk['chr'] == chromosome] for chunk in snps_df_iter])
    snps_df.drop('chr', axis=1)
    snps_df.set_index('pos', inplace=True)

    snp_counts_df = snp_counts_df.merge(snps_df, left_index=True, right_index=True)

    # Create a list of heterozygous snps to search for in each tumour
    het_df = snp_counts_df[snp_counts_df['AB'] == 1]
    het_df.reset_index(inplace=True)
    het_df['chr'] = chromosome
    het_df.to_csv(hets_filename, sep='\t', cols=['chr', 'pos', 'ref', 'alt'], index=False, header=False)

    # Create genotype file required by impute2
    temp_gen_filename = os.path.join(haps_temp_directory, 'snps.gen')
    snp_counts_df.reset_index(inplace=True)
    snp_counts_df['chr'] = chromosome
    snp_counts_df['chr_pos'] = snp_counts_df['chr'].astype(str) + ':' + snp_counts_df['pos'].astype(str)
    snp_counts_df.to_csv(temp_gen_filename, sep=' ', cols=['chr', 'chr_pos', 'pos', 'ref', 'alt', 'AA', 'AB', 'BB'], index=False, header=False)

    # Create single sample file required by impute2
    temp_sample_filename = os.path.join(haps_temp_directory, 'snps.sample')
    with open(temp_sample_filename, 'w') as temp_sample_file:
        temp_sample_file.write('ID_1 ID_2 missing sex\n0 0 0 0\nUNR1 UNR1 0 2\n')

    # Run shapeit to create phased haplotype graph
    hgraph_filename = os.path.join(haps_temp_directory, 'phased.hgraph')
    hgraph_logs_prefix = hgraph_filename + '.log'
    chr_x_flag = ''
    if chromosome == 'X':
        chr_x_flag = '--chrX'
    pypeliner.commandline.execute('shapeit', '-M', genetic_map_filename, '-R', hap_filename, legend_filename, config['sample_filename'],
                                  '-G', temp_gen_filename, temp_sample_filename, '--output-graph', hgraph_filename, chr_x_flag,
                                  '--no-mcmc', '-L', hgraph_logs_prefix)

    # Run shapeit to sample from phased haplotype graph
    sample_template = os.path.join(haps_temp_directory, 'sampled.{0}')
    averaged_changepoints = None
    for s in range(int(config['shapeit_num_samples'])):
        sample_prefix = sample_template.format(s)
        sample_log_filename = sample_prefix + '.log'
        sample_haps_filename = sample_prefix + '.haps'
        sample_sample_filename = sample_prefix + '.sample'
        pypeliner.commandline.execute('shapeit', '-convert', '--input-graph', hgraph_filename, '--output-sample', 
                                      sample_prefix, '--seed', str(s), '-L', sample_log_filename)
        sample_haps = pd.read_csv(sample_haps_filename, sep=' ', header=None, 
                                  names=['id', 'id2', 'pos', 'ref', 'alt', 'allele1', 'allele2'],
                                  usecols=['pos', 'allele1', 'allele2'])
        sample_haps = sample_haps[sample_haps['allele1'] != sample_haps['allele2']]
        sample_haps['allele'] = sample_haps['allele1']
        sample_haps = sample_haps.drop(['allele1', 'allele2'], axis=1)
        sample_haps.set_index('pos', inplace=True)
        sample_changepoints = sample_haps['allele'].diff().abs().astype(float).fillna(0.0)
        if averaged_changepoints is None:
            averaged_changepoints = sample_changepoints
        else:
            averaged_changepoints += sample_changepoints
        os.remove(sample_log_filename)
        os.remove(sample_haps_filename)
        os.remove(sample_sample_filename)
    averaged_changepoints /= float(config['shapeit_num_samples'])
    last_sample_haps = sample_haps

    # Identify changepoints recurrent across samples
    changepoint_confidence = np.maximum(averaged_changepoints, 1.0 - averaged_changepoints)

    # Create a list of labels for haplotypes between recurrent changepoints
    current_hap_label = 0
    hap_label = list()
    for x in changepoint_confidence:
        if x < float(config['shapeit_confidence_threshold']):
            current_hap_label += 1
        hap_label.append(current_hap_label)

    # Create the list of haplotypes
    haps = last_sample_haps
    haps['changepoint_confidence'] = changepoint_confidence
    haps['hap_label'] = hap_label

    haps.reset_index(inplace=True)

    haps['allele_id'] = 0

    haps_allele2 = haps.copy()
    haps_allele2['allele_id'] = 1
    haps_allele2['allele'] = 1 - haps_allele2['allele']

    haps = pd.concat([haps, haps_allele2], ignore_index=True)
    haps.sort(['pos', 'allele_id'], inplace=True)

    haps.set_index(['hap_label', 'allele_id'], inplace=True)
    hap_label_counter = itertools.count()
    haps['allele_label'] = haps.groupby(level=[0, 1]).apply(lambda a: next(hap_label_counter))
    haps.reset_index(inplace=True)

    haps = haps[['pos', 'allele', 'changepoint_confidence', 'hap_label', 'allele_id', 'allele_label']]

    haps.to_csv(haps_filename, sep='\t', header=True, index=False)


def merge_files(output_filename, *input_filenames):
    with open(output_filename, 'w') as output_file:
        for input_filename in input_filenames:
            with open(input_filename, 'r') as input_file:
                shutil.copyfileobj(input_file, output_file)


def merge_tables(output_filename, *input_filenames):
    output_table = list()
    for input_filename in input_filenames:
        output_table.append(pd.read_csv(input_filename, sep='\t'))
    output_table = pd.concat(output_table, ignore_index=True)
    output_table.to_csv(output_filename, sep='\t', index=False)


def merge_tars(output_filename, *input_filenames):
    with tarfile.open(output_filename, 'w') as output_tar:
        for input_filename in input_filenames:
            with tarfile.open(input_filename, 'r') as in_tar:
                for tarinfo in in_tar:
                    output_tar.addfile(tarinfo, in_tar.extractfile(tarinfo))


def read_chromosome_lengths(genome_fai_filename):
    chromosome_lengths = dict()
    with open(genome_fai_filename, 'r') as genome_fai_file:
        for row in csv.reader(genome_fai_file, delimiter='\t'):
            chromosome = row[0]
            length = int(row[1])
            if chromosome.startswith('GL'):
                continue
            if chromosome == 'MT':
                continue
            chromosome_lengths[chromosome] = length
    return chromosome_lengths


def filled_density(ax, data, c, xmin, xmax, cov):
    density = scipy.stats.gaussian_kde(data)
    density.covariance_factor = lambda : cov
    density._compute_covariance()
    xs = [xmin] + list(np.linspace(xmin, xmax, 2000)) + [xmax]
    ys = density(xs)
    ys[0] = 0.0
    ys[-1] = 0.0
    ax.plot(xs, ys, color=c)
    return ax.fill(xs, ys, color=c, alpha=0.5)


def savefig_tar(tar, fig, filename):
    plot_buffer = io.BytesIO()
    fig.savefig(plot_buffer, format='pdf')
    info = tarfile.TarInfo(name=filename)
    info.size = plot_buffer.tell()
    plot_buffer.seek(0)
    tar.addfile(tarinfo=info, fileobj=plot_buffer)


def phase_intervals(*args):

    input_alleles_filenames = args[:len(args)/2]
    output_alleles_filenames = args[len(args)/2:]

    allele_phases = list()
    allele_diffs = list()

    for idx, input_alleles_filename in enumerate(input_alleles_filenames):

        allele_data = pd.read_csv(input_alleles_filename, sep='\t', header=None, names=['interval_id', 'hap_label', 'allele_id', 'readcount'])
        
        # Allele readcount table
        allele_data = allele_data.set_index(['interval_id', 'hap_label', 'allele_id'])['readcount'].unstack().fillna(0.0)
        allele_data = allele_data.astype(float)
        
        # Create major allele call
        allele_phase = allele_data.apply(np.argmax, axis=1)
        allele_phase.name = 'major_allele_id'
        allele_phase = allele_phase.reset_index()
        allele_phase['library_idx'] = idx
        allele_phases.append(allele_phase)

        # Calculate major minor allele read counts, and diff between them
        allele_data['major_readcount'] = allele_data.apply(np.max, axis=1)
        allele_data['minor_readcount'] = allele_data.apply(np.min, axis=1)
        allele_data['diff_readcount'] = allele_data['major_readcount'] - allele_data['minor_readcount']
        allele_data['total_readcount'] = allele_data['major_readcount'] + allele_data['minor_readcount']

        # Calculate normalized major and minor read counts difference per interval
        allele_diff = allele_data.groupby(level=[0])[['diff_readcount', 'total_readcount']].sum()
        allele_diff['norm_diff_readcount'] = allele_diff['diff_readcount'] / allele_diff['total_readcount']
        allele_diff = allele_diff[['norm_diff_readcount']]

        # Add to table for all librarys
        allele_diff.reset_index(inplace=True)
        allele_diff['library_idx'] = idx
        allele_diffs.append(allele_diff)

    allele_phases = pd.concat(allele_phases, ignore_index=True)
    allele_diffs = pd.concat(allele_diffs, ignore_index=True)

    def select_largest_diff(df):
        largest_idx = np.argmax(df['norm_diff_readcount'].values)
        return df['library_idx'].values[largest_idx]

    # For each interval, select the library with the largest difference between major and minor
    interval_library = allele_diffs.set_index('interval_id').groupby(level=0).apply(select_largest_diff)
    interval_library.name = 'library_idx'
    interval_library = interval_library.reset_index()

    # For each haplotype block in each interval, take the major allele call of the library
    # with the largest major minor difference and call it allele 'a'
    allele_phases = allele_phases.merge(interval_library, left_on=['interval_id', 'library_idx'], right_on=['interval_id', 'library_idx'], how='right')
    allele_phases = allele_phases[['interval_id', 'hap_label', 'major_allele_id']].rename(columns={'major_allele_id': 'allele_a_id'})

    for idx, (input_alleles_filename, output_allele_filename) in enumerate(zip(input_alleles_filenames, output_alleles_filenames)):

        allele_data = pd.read_csv(input_alleles_filename, sep='\t', header=None, names=['interval_id', 'hap_label', 'allele_id', 'readcount'])
        
        # Add a boolean column denoting which allele is allele 'a'
        allele_data = allele_data.merge(allele_phases, left_on=['interval_id', 'hap_label'], right_on=['interval_id', 'hap_label'])
        allele_data['is_allele_a'] = (allele_data['allele_id'] == allele_data['allele_a_id']) * 1
        allele_data = allele_data[['interval_id', 'hap_label', 'allele_id', 'readcount', 'is_allele_a']]

        allele_data.to_csv(output_allele_filename, sep='\t', header=False, index=False)


def gc_lowess(gc_samples_filename, gc_dist_filename, plot_filename, gc_resolution=100):

    gc_samples = pd.read_csv(gc_samples_filename, sep='\t', names=['chromosome', 'position', 'gc', 'count'])

    gc_samples['gc_bin'] = np.round(gc_samples['gc'] * gc_resolution)

    gc_binned = gc_samples.groupby('gc_bin')['count'] \
                          .agg({'sum':np.sum, 'len':len, 'mean':np.mean}) \
                          .reindex(xrange(gc_resolution+1)) \
                          .fillna(0) \
                          .reset_index() \
                          .rename(columns={'index':'gc_bin'}) \
                          .astype(float)

    gc_binned['smoothed'] = sm.nonparametric.lowess(gc_binned['mean'].values, gc_binned['gc_bin'].values, frac=0.1).T[1]

    rescale = 1. / gc_binned['smoothed'].max()

    gc_binned['mean'] = gc_binned['mean'] * rescale
    gc_binned['smoothed'] = gc_binned['smoothed'] * rescale

    fig = plt.figure(figsize=(4,4))

    plt.scatter(gc_binned['gc_bin'].values, gc_binned['mean'].values, c='k', s=4)
    plt.plot(gc_binned['gc_bin'].values, gc_binned['smoothed'].values, c='r')

    plt.xlabel('gc %')
    plt.ylabel('density')
    plt.xlim((-0.5, 100.5))
    plt.ylim((-0.01, gc_binned['mean'].max() * 1.1))

    plt.tight_layout()

    fig.savefig(plot_filename, format='pdf', bbox_inches='tight')

    gc_binned[['smoothed']].to_csv(gc_dist_filename, sep='\t', index=False, header=False)


def solve_and_plot(library_id, intervals_filename, alleles_filename, stats_filename, pred_filename, plots_tar_filename, plots_prefix):

    with tarfile.open(plots_tar_filename, 'w') as plots_tar:

        high_confidence_length = 100000.0

        interval_data = pd.read_csv(intervals_filename, sep='\t', header=None, converters={'id':str, 'chromosome1':str, 'chromosome2':str},
                                    names=['id', 'chromosome1', 'position1', 'strand1', 'chromosome2', 'position2', 'strand2', 'readcount', 'length'])

        allele_data = pd.read_csv(alleles_filename, sep='\t', header=None, names=['interval_id', 'hap_label', 'allele_id', 'readcount', 'is_allele_a'])

        # Calculate allele a/b readcounts
        allele_data = allele_data.set_index(['interval_id', 'hap_label', 'is_allele_a'])['readcount'].unstack().fillna(0.0)
        allele_data = allele_data.astype(int)
        allele_data = allele_data.rename(columns={0:'allele_b_readcount', 1:'allele_a_readcount'})

        # Calculate expected read depth based on total read count and total length
        allele_data.reset_index(inplace=True)
        allele_data.set_index('interval_id', inplace=True)
        interval_data.set_index('id', inplace=True)
        allele_data['expected_depth'] = (interval_data['readcount'] / interval_data['length'])
        allele_data.reset_index(inplace=True)

        # Estimate effective haplotype block length based on expected depth and total allele readcounts
        allele_data.set_index(['interval_id', 'hap_label'], inplace=True)
        allele_data['hap_length'] = (allele_data['allele_b_readcount'] + allele_data['allele_a_readcount']) / allele_data['expected_depth']
        allele_data = allele_data.replace([np.inf, -np.inf], np.nan).dropna()

        # Merge haplotype blocks contained within the same interval
        allele_data = allele_data.groupby(level=[0])[['hap_length', 'allele_a_readcount', 'allele_b_readcount']].sum()

        # Calculate major and minor readcounts, and relationship to allele a/b
        allele_data['major_readcount'] = allele_data[['allele_a_readcount', 'allele_b_readcount']].apply(max, axis=1)
        allele_data['minor_readcount'] = allele_data[['allele_a_readcount', 'allele_b_readcount']].apply(min, axis=1)
        allele_data['major_is_allele_a'] = (allele_data['major_readcount'] == allele_data['allele_a_readcount']) * 1

        # Merge allele data with interval data
        interval_data = interval_data.merge(allele_data, left_index=True, right_index=True)

        # Calculate coverages
        interval_data['minor_cov'] = interval_data['minor_readcount'] / interval_data['hap_length']
        interval_data['major_cov'] = interval_data['major_readcount'] / interval_data['hap_length']
        interval_data['total_cov'] = interval_data['readcount'] / interval_data['length']
        interval_data = interval_data.fillna(0.0)

        # Classify as high or low confidence depending on length and haplotype length
        interval_data['high_conf'] = 0
        interval_data.loc[(interval_data['length'] > high_confidence_length) & (interval_data['hap_length'] > 0), 'high_conf'] = 1

        # Split off the low confidence intervals
        interval_data_low_conf = interval_data.loc[interval_data['high_conf'] == 0]
        interval_data = interval_data.loc[interval_data['high_conf'] == 1]

        # Sort by chromosome / position, required for variance estimates
        interval_data = interval_data.sort(['chromosome1', 'position1'])

        # Obtain major and minor allele depth for high confidence intervals
        minor = interval_data['minor_cov'].values
        major = interval_data['major_cov'].values
        total = interval_data['total_cov'].values

        # Other parameters from high confidence intervals
        chrs = interval_data['chromosome1'].values
        starts = interval_data['position1'].values
        ends = interval_data['position2'].values
        lengths = interval_data['length'].values
        major_is_allele_a = interval_data['major_is_allele_a'].values

        # Calculate difference between adjacent intervals
        minor_sq_diffs = np.square(minor[1:] - minor[:-1])
        major_sq_diffs = np.square(major[1:] - major[:-1])
        length_sums = lengths[1:] + lengths[:-1]

        # Resample squared differences by sum of lengths of adjacent intervals
        samples_length_sums = np.array(list(itertools.chain(*[itertools.repeat(idx, cnt) for idx, cnt in enumerate(np.random.multinomial(100000, length_sums / length_sums.sum()))])))
        minor_sq_diff_samples = minor_sq_diffs[samples_length_sums]
        major_sq_diff_samples = major_sq_diffs[samples_length_sums]

        # Estimate variance as mean squared difference between adjacent regions
        minor_variance_estimate = np.median(minor_sq_diff_samples)
        major_variance_estimate = np.median(major_sq_diff_samples)

        # Resample major and minor allele depths weighted by lengths
        samples_length = np.array(list(itertools.chain(*[itertools.repeat(idx, cnt) for idx, cnt in enumerate(np.random.multinomial(10000, lengths / lengths.sum()))])))
        minor_samples = minor[samples_length]
        major_samples = major[samples_length]
        total_samples = total[samples_length]

        # Maximum depth for plots
        depth_max = np.percentile(major_samples + minor_samples, 95.0)

        def gmm_bic(n_comp, X):
            kmm = sklearn.cluster.KMeans(n_clusters=n_comp)
            kmm.fit(X.reshape((X.size, 1)))
            mm = sklearn.mixture.GMM(n_components=n_comp, covariance_type='tied', n_iter=int(1e6), thresh=1e-9, init_params='w', params='mw')
            mm.means_ = kmm.cluster_centers_
            mm.covars_ = np.array([minor_variance_estimate]).reshape((1,1))
            mm.fit(X)
            assert mm.converged_
            return mm.bic(X), mm.means_, mm.weights_

        # Fit a GMM to normal read depths using BIC
        bic, means, weights = sorted(gmm_bic(a, minor_samples) for a in range(1, 12+1))[0]
        means = means.reshape((means.size,))
        weights = weights.reshape((weights.size,))

        # Use modes in minor depth as starting points for normal contamination and haploid coverage
        weight_filter = weights >= 0.05
        normal_mean = min(means[weight_filter])
        tumour_means = means[weight_filter & (means != normal_mean)] - normal_mean
        unique_tumour_means = np.unique(np.around(tumour_means, 6))
        unique_tumour_means = unique_tumour_means[unique_tumour_means > 0.0]
        haploid_tumour_candidates = np.concatenate((unique_tumour_means, 0.5*unique_tumour_means))

        # Fallback if no significant second mode
        if len(haploid_tumour_candidates) == 0:
            haploid_tumour_candidates = np.array([normal_mean * 0.5])

        # Plot raw alleles density with GMM modes
        fig = plt.figure(figsize=(16,8))
        filled_density(plt.gca(), minor_samples, 'blue', 0.0, depth_max, 0.01)
        filled_density(plt.gca(), major_samples, 'red', 0.0, depth_max, 0.01)
        filled_density(plt.gca(), total_samples, 'grey', 0.0, depth_max, 0.01)
        x = np.linspace(-5.0*np.sqrt(minor_variance_estimate), 5.0*np.sqrt(minor_variance_estimate), 1000.0)
        plt.plot(x, scipy.stats.norm.pdf(x, scale=np.sqrt(minor_variance_estimate)), 'blue')
        x = np.linspace(-5.0*np.sqrt(major_variance_estimate), 5.0*np.sqrt(major_variance_estimate), 1000.0)
        plt.plot(x, scipy.stats.norm.pdf(x, scale=np.sqrt(major_variance_estimate)), 'red')
        ylim = plt.ylim()
        for mean in means:
            if mean == normal_mean:
                plt.plot([mean, mean], [0, 1e16], 'g', lw=2)
            else:
                plt.plot([mean, mean], [0, 1e16], 'g--', lw=2)
        plt.xlim((0.0, depth_max))
        plt.ylim(ylim)
        savefig_tar(plots_tar, fig, '{0}_{1}_raw_allele_density.pdf'.format(plots_prefix, library_id))
        plt.clf()

        model_probs = list()
        haploid_normals = list()
        haploid_tumours = list()
        subclone_freqs = list()
        avg_zs = list()
        preds = list()

        # Iterate over potential haploid coverage candidates
        for candidate_mean in haploid_tumour_candidates:

            # Estimate normal contamination and haploid tumour coverage
            model_prob, haploid_normal, haploid_tumour, subclone_freq, avg_z, pred = subclonal_sampling.estimate_model_prob(major, minor, lengths, normal_mean, candidate_mean, major_variance_estimate, minor_variance_estimate)

            model_probs.append(model_prob)
            haploid_normals.append(haploid_normal)
            haploid_tumours.append(haploid_tumour)
            subclone_freqs.append(subclone_freq)
            avg_zs.append(avg_z)
            preds.append(pred)

        def close_enough(x, y, max_diff=0.01):
            return abs(x - y) / (0.5 * (x + y)) < max_diff

        # Create list of unique predictions
        filtered = list(range(len(model_probs)))
        for idx1 in xrange(len(model_probs)):
            for idx2 in xrange(len(model_probs)):
                if close_enough(haploid_normals[idx1], haploid_normals[idx2]) and close_enough(haploid_tumours[idx1], haploid_tumours[idx2]) and model_probs[idx1] > model_probs[idx2]:
                    try:
                        filtered.remove(idx2)
                    except:
                        pass

        # Filter non unique predictions
        model_probs = np.array(model_probs)[filtered]
        haploid_normals = np.array(haploid_normals)[filtered]
        haploid_tumours = np.array(haploid_tumours)[filtered]
        subclone_freqs = np.array(subclone_freqs)[filtered]
        avg_zs = np.array(avg_zs)[filtered]
        preds = np.array(preds)[filtered]

        stats_tables = list()
        preds_tables = list()

        for candidate_idx, (model_prob, haploid_normal, haploid_tumour, subclone_freq, avg_z, pred) in enumerate(zip(model_probs, haploid_normals, haploid_tumours, subclone_freqs, avg_zs, preds)):

            # Shift modes by estimated normal contamination and scale by haploid coverage
            tumour_minor_samples = (minor_samples - haploid_normal) / haploid_tumour
            tumour_major_samples = (major_samples - haploid_normal) / haploid_tumour
            tumour_depth_max = (depth_max - haploid_normal) / haploid_tumour
            tumour_minor = (minor - haploid_normal) / haploid_tumour
            tumour_major = (major - haploid_normal) / haploid_tumour

            # Extract predictions
            pred_major = pred[0]
            pred_minor = pred[1]
            pred_major_sub = pred[2]
            pred_minor_sub = pred[3]

            # Plot adjusted alleles density with copy number calls
            fig = plt.figure(figsize=(16,8))
            filled_density(plt.gca(), tumour_minor_samples, 'blue', -0.5, tumour_depth_max, 0.01)
            filled_density(plt.gca(), tumour_major_samples, 'red', -0.5, tumour_depth_max, 0.01)
            filled_density(plt.gca(), tumour_minor_samples + tumour_major_samples, 'grey', -0.5, tumour_depth_max, 0.01)
            ylim = plt.ylim()
            for i in range(20):
                if i == 0:
                    plt.plot([i, i], [0, 1e16], 'g', lw=2)
                else:
                    plt.plot([i, i], [0, 1e16], 'g--', lw=2)
            plt.xlim((-0.5, tumour_depth_max))
            plt.ylim(ylim)
            savefig_tar(plots_tar, fig, '{0}_{1}_{2}_allele_density.pdf'.format(plots_prefix, library_id, candidate_idx))
            plt.clf()

            # Plot major vs minor tumour copies annotated by chromosome
            chromosomes = [str(a) for a in range(1, 23)] + ['X']
            color_set = plt.get_cmap('Set1')
            color_set = [color_set(float(i)/len(chromosomes)) for i in range(len(chromosomes))]
            chromosome_color = lambda c: color_set[chromosomes.index(c)]
            cs = [chromosome_color(c) for c in chrs]
            fig = plt.figure(figsize=(20,16))
            plt.scatter(tumour_major, tumour_minor, s=lengths/20000.0, facecolor=cs, edgecolor=cs, linewidth=0.0)
            plt.grid(True)
            plt.legend([plt.Circle((0, 0), radius=1, color=chromosome_color(c)) for c in chromosomes], chromosomes, loc=2)
            plt.xlim((-0.5, tumour_depth_max))
            plt.ylim((-0.5, 0.8*tumour_depth_max))
            plt.xlabel('major copy number')
            plt.ylabel('minor copy number')
            savefig_tar(plots_tar, fig, '{0}_{1}_{2}_major_minor_chromosome.pdf'.format(plots_prefix, library_id, candidate_idx))
            plt.clf()

            # Plot major vs minor tumour copies annotated by subclonality
            color_set = plt.get_cmap('copper')
            cs = [color_set(z*2.) for z in avg_z]
            fig = plt.figure(figsize=(20,16))
            plt.scatter(tumour_major, tumour_minor, s=lengths/20000.0, facecolor=cs, edgecolor=cs, linewidth=0.0)
            plt.grid(True)
            plt.xlim((-0.5, tumour_depth_max))
            plt.ylim((-0.5, 0.8*tumour_depth_max))
            plt.xlabel('major copy number')
            plt.ylabel('minor copy number')
            savefig_tar(plots_tar, fig, '{0}_{1}_{2}_major_minor_subclonality.pdf'.format(plots_prefix, library_id, candidate_idx))
            plt.clf()

            # Calculate tumour genome size (assuming human normal genome size 3e9)
            minor_size = np.sum(tumour_minor * lengths) * 3e9 / np.sum(lengths)
            major_size = np.sum(tumour_major * lengths) * 3e9 / np.sum(lengths)

            # Create stats table
            stats = pd.DataFrame({'library_id':[library_id]})
            stats['candidate_id'] = candidate_idx
            stats['likelihood'] = model_prob
            stats['haploid_normal'] = haploid_normal
            stats['haploid_tumour'] = haploid_tumour
            stats['minor_size'] = minor_size
            stats['major_size'] = major_size
            stats['minor_variance'] = minor_variance_estimate
            stats['major_variance'] = major_variance_estimate
            stats['tumour_cell_proportion'] = stats['haploid_tumour'] / (stats['haploid_tumour'] + stats['haploid_normal'])
            stats['normal_cell_proportion'] = 1.0 - stats['tumour_cell_proportion']
            stats['subclone_frequency'] = subclone_freq
            stats_tables.append(stats)

            # Create copy number predictions
            preds = pd.DataFrame({'chr':chrs, 'start':starts, 'end':ends, 'length':lengths, 'major_raw':tumour_major, 'minor_raw':tumour_minor, 'major':pred_major, 'minor':pred_minor, 'major_sub':pred_major_sub, 'minor_sub':pred_minor_sub, 'subclonal':avg_z, 'major_is_allele_a':major_is_allele_a})
            preds['high_conf'] = 1

            # Predictions for low confidence intervals
            low_conf_preds = list()
            for idx, row in interval_data_low_conf.iterrows():
                low_conf_preds.append(subclonal_sampling.max_major_minor_posterior(row['major_cov'], row['minor_cov'], haploid_normal, haploid_tumour, 0.0, major_variance_estimate, minor_variance_estimate))
            low_conf_preds = pd.DataFrame(low_conf_preds, columns=['major', 'minor', 'major_sub', 'minor_sub'])
            low_conf_preds['chr'] = interval_data_low_conf['chromosome1'].values
            low_conf_preds['start'] = interval_data_low_conf['position1'].values
            low_conf_preds['end'] = interval_data_low_conf['position2'].values
            low_conf_preds['length'] = interval_data_low_conf['length'].values
            low_conf_preds['major_raw'] = ((interval_data_low_conf['major_cov'] - haploid_normal) / haploid_tumour).values
            low_conf_preds['minor_raw'] = ((interval_data_low_conf['minor_cov'] - haploid_normal) / haploid_tumour).values
            low_conf_preds['major_is_allele_a'] = interval_data_low_conf['major_is_allele_a'].values
            low_conf_preds['high_conf'] = 0

            # Combine low and high confidence
            preds = pd.concat([preds, low_conf_preds], axis=0, ignore_index=True).sort(['chr', 'start'])

            # Add lib id and candidate id and append to list of tables
            preds['library_id'] = library_id
            preds['candidate_id'] = candidate_idx
            preds = preds[['library_id', 'candidate_id', 'chr', 'start', 'end', 'length', 'major', 'minor', 'major_sub', 'minor_sub', 'subclonal', 'major_raw', 'minor_raw', 'high_conf', 'major_is_allele_a']]
            preds_tables.append(preds)

        pd.concat(stats_tables, ignore_index=True).to_csv(stats_filename, sep='\t', na_rep='NA', index=False, header=True)
        pd.concat(preds_tables, ignore_index=True).to_csv(pred_filename, sep='\t', na_rep='NA', index=False, header=True)

