import csv
import os

import pypeliner

import destruct.defaultconfig


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
        with open(sentinal_filename, 'w'):
            pass

def create_ref_data(config, ref_data_dir):
    config = destruct.defaultconfig.get_config(ref_data_dir, config)

    try:
        os.makedirs(ref_data_dir)
    except OSError:
        pass

    auto_sentinal = AutoSentinal(ref_data_dir + '/sentinal.')

    temp_directory = os.path.join(ref_data_dir, 'tmp')

    try:
        os.makedirs(temp_directory)
    except OSError:
        pass

    def wget_genome_fasta():
        with open(config['genome_fasta'], 'w') as genome_file:
            for assembly in config['ensembl_assemblies']:
                assembly_url = config['ensembl_assembly_url'].format(assembly)
                assembly_fasta = os.path.join(temp_directory, 'dna.assembly.{0}.fa'.format(assembly))
                if not os.path.exists(assembly_fasta):
                    wget_gunzip(assembly_url, assembly_fasta)
                with open(assembly_fasta, 'r') as assembly_file:
                    for line in assembly_file:
                        if line[0] == '>':
                            line = line.split()[0] + '\n'
                        genome_file.write(line)
    auto_sentinal.run(wget_genome_fasta)

    def wget_gtf():
        wget_gunzip(config['ensembl_gtf_url'], config['gtf_filename'])
    auto_sentinal.run(wget_gtf)

    def wget_dgv():
        wget(config['dgv_url'], config['dgv_filename'])
    auto_sentinal.run(wget_dgv)

    def wget_repeats():
        repeat_filename = os.path.join(temp_directory, 'repeats.txt')
        wget_gunzip(config['rmsk_url'], repeat_filename)
        with open(config['chromosome_map'], 'r') as chr_map_file:
            chr_map = dict((a.split() for a in chr_map_file))
        with open(repeat_filename, 'r') as repeat_file, open(config['repeat_regions'], 'w') as repeat_regions_file, open(config['satellite_regions'], 'w') as satellite_regions_file:
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
        pypeliner.commandline.execute('bowtie-build', config['genome_fasta'], config['genome_fasta'])
    auto_sentinal.run(bowtie_build)

    def samtools_faidx():
        pypeliner.commandline.execute('samtools', 'faidx', config['genome_fasta'])
    auto_sentinal.run(samtools_faidx)


