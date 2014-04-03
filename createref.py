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
import gzip
from collections import *

import pypeliner

def wget_gunzip(url, filename):
    temp_filename = filename + '.tmp'
    pypeliner.commandline.execute('wget', url, '-c', '-O', temp_filename + '.gz')
    pypeliner.commandline.execute('gunzip', temp_filename + '.gz')
    os.rename(temp_filename, filename)

def wget(url, filename):
    temp_filename = filename + '.tmp'
    pypeliner.commandline.execute('wget', url, '-c', '-O', temp_filename)
    os.rename(temp_filename, filename)

class AutoSentinal(object):
    def __init__(self, sentinal_prefix):
        self.sentinal_prefix = sentinal_prefix
    def run(self, func):
        sentinal_filename = self.sentinal_prefix + func.__name__
        if os.path.exists(sentinal_filename):
            return
        func()
        with open(sentinal_filename, 'w') as sentinal_file:
            pass

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('config', help='Configuration filename')
    parser.add_argument('-d', '--demix', action='store_true', help='Download additional demix data')
    args = parser.parse_args()

    args.config = os.path.abspath(args.config)

    config = ConfigParser.ConfigParser()
    config.read(args.config)

    class ConfigWrapper(object):
        def __init__(self, config, section='main'):
            self._config = config
            self._section = section
        def __getattr__(self, name):
            if name.endswith('_tool'):
                return self._config.get(self._section, 'tools_directory') + '/' + name[:-5]
            if name.endswith('_rscript'):
                return self._config.get(self._section, 'rscripts_directory') + '/' + name[:-8] + '.R'
            return self._config.get(self._section, name)

    cfg = ConfigWrapper(config)

    try:
        os.makedirs(cfg.dataset_directory)
    except OSError:
        pass

    auto_sentinal = AutoSentinal(cfg.dataset_directory + '/sentinal.')


    def wget_genome_fasta():
        assembly_prefix = cfg.dataset_directory + '/dna.assembly.'
        assembly_fastas = dict()
        for assembly in cfg.ensembl_assemblies.rstrip().split(','):
            assembly_fastas[assembly] = assembly_prefix + assembly + '.fa'
            if os.path.exists(assembly_fastas[assembly]):
                continue
            assembly_url = 'ftp://ftp.ensembl.org/pub/release-{0}/fasta/homo_sapiens/dna/Homo_sapiens.{1}.{2}.dna.{3}.fa.gz'.format(cfg.ensembl_version, cfg.ensembl_genome_version, cfg.ensembl_version, assembly)
            wget_gunzip(assembly_url, assembly_fastas[assembly])
            with open(cfg.genome_fasta, 'w') as genome_file:
                for assembly, assembly_fasta in assembly_fastas.iteritems():
                    with open(assembly_fasta, 'r') as assembly_file:
                        for line in assembly_file:
                            if line[0] == '>':
                                line = line.split()[0] + '\n'
                            genome_file.write(line)
    auto_sentinal.run(wget_genome_fasta)

    def wget_gtf():
        gtf_url = 'ftp://ftp.ensembl.org/pub/release-{0}/gtf/homo_sapiens/Homo_sapiens.{1}.{2}.gtf.gz'.format(cfg.ensembl_version, cfg.ensembl_genome_version, cfg.ensembl_version)
        wget_gunzip(gtf_url, cfg.gtf_filename)
    auto_sentinal.run(wget_gtf)

    def wget_dgv():
        if cfg.dgv_url != '':
            dgv_url = cfg.dgv_url
        else:
            dgv_url = 'http://dgv.tcag.ca/dgv/docs/{0}_variants_{1}.txt'.format(cfg.dgv_genome_version, cfg.dgv_version)
        wget(dgv_url, cfg.dgv_filename)
    auto_sentinal.run(wget_dgv)

    def wget_repeats():
        repeat_filename = cfg.dataset_directory + '/repeats.txt'
        if cfg.ucsc_genome_version == 'hg18':
            rmsk_url = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/*_rmsk.txt.gz'
            wget_gunzip(rmsk_url, repeat_filename)
        else:
            rmsk_url = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/{0}/database/rmsk.txt.gz'.format(cfg.ucsc_genome_version)
            wget_gunzip(rmsk_url, repeat_filename)
        with open(cfg.chromosome_map, 'r') as chr_map_file:
            chr_map = dict((a.split() for a in chr_map_file))
        with open(repeat_filename, 'r') as repeat_file, open(cfg.repeat_regions, 'w') as repeat_regions_file, open(cfg.satellite_regions, 'w') as satellite_regions_file:
            for row in csv.reader(repeat_file, delimiter='\t'):
                if row[5] not in chr_map:
                    continue
                chr = chr_map[row[5]]
                start = str(int(row[6]) + 1)
                end = row[7]
                type = row[11]
                repeat_regions_file.write('\t'.join([chr, start, end]) + '\n')
                if type == 'Satellite':
                    satellite_regions_file.write('\t'.join([chr, start, end]) + '\n')
    auto_sentinal.run(wget_repeats)

    def bowtie_build():
        pypeliner.commandline.execute(cfg.bowtie_build_bin, cfg.genome_fasta, cfg.genome_fasta)
    auto_sentinal.run(bowtie_build)

    def bowtie2_build():
        pypeliner.commandline.execute(cfg.bowtie2_build_bin, cfg.genome_fasta, cfg.genome_fasta)
    auto_sentinal.run(bowtie2_build)

    def mrsfast_index():
        pypeliner.commandline.execute(cfg.mrsfast_bin, '--index', cfg.genome_fasta)
    auto_sentinal.run(mrsfast_index)

    def samtools_faidx():
        pypeliner.commandline.execute(cfg.samtools_bin, 'faidx', cfg.genome_fasta)
    auto_sentinal.run(samtools_faidx)


    if args.demix:

        def wget_thousand_genomes():
            tar_filename = os.path.join(cfg.dataset_directory, 'thousand_genomes_download.tar.gz')
            wget(cfg.thousand_genomes_impute_url, tar_filename)
            pypeliner.commandline.execute('tar', '-C', cfg.dataset_directory, '-xzvf', tar_filename)
            os.remove(tar_filename)
        auto_sentinal.run(wget_thousand_genomes)

        def create_snp_positions():
            with open(cfg.snp_positions, 'w') as snp_positions_file:
                for chromosome in cfg.chromosomes.split():
                    phased_chromosome = chromosome
                    if chromosome == 'X':
                        phased_chromosome = cfg.phased_chromosome_x
                    legend_filename = cfg.legend_template.format(phased_chromosome)
                    with gzip.open(legend_filename, 'r') as legend_file:
                        for line in legend_file:
                            if line.startswith('id'):
                                continue
                            row = line.split()
                            rs_id = row[0]
                            position = row[1]
                            a0 = row[2]
                            a1 = row[3]
                            if len(a0) != 1 or len(a1) != 1:
                                continue
                            snp_positions_file.write('\t'.join([chromosome, position, a0, a1]) + '\n')
        auto_sentinal.run(create_snp_positions)

