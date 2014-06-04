import argparse
import ConfigParser
import random
import string
import os
import subprocess
import sys


def read_sequences(fasta):
    id = None
    sequences = []
    for line in fasta:
        line = line.rstrip()
        if len(line) == 0:
            continue
        if line[0] == '>':
            if id is not None:
                yield (id, ''.join(sequences))
            id = line[1:]
            sequences = []
        else:
            sequences.append(line)
    if id is not None:
        yield (id, ''.join(sequences))

def reverse_complement(sequence):
    return sequence[::-1].translate(string.maketrans('ACTGactg','TGACtgac'))

def random_chromosome_position(genome, dist_to_end):
    while True:
        position = random.randint(1, sum([len(chr) for chr in genome.itervalues()]))
        for id, chr in genome.iteritems():
            if position <= len(chr):
                if position <= dist_to_end or len(chr) - position < dist_to_end:
                    break
                return id, random.choice(['+', '-']), position
            position -= len(chr)

def retrieve_upstream_sequence(genome, chr, str, pos, length):
    pos = pos - 1
    if str == '+':
        return genome[chr][pos-length+1:pos+1]
    else:
        return genome[chr][pos:pos+length]

def retrieve_downstream_sequence(genome, chr, str, pos, length):
    pos = pos - 1
    if str == '+':
        return genome[chr][pos+1:pos+length+1]
    else:
        return genome[chr][pos-length:pos]

def create_distinct(sequence1, sequence2, length):
    sequence1 = sequence1.upper()
    sequence2 = sequence2.upper()
    distinct = ''
    for idx, nts in enumerate(zip(sequence1, sequence2)):
        while len(distinct) < length:
            nt = random.choice(['A', 'C', 'T', 'G'])
            if nt not in nts:
                distinct = distinct + nt
                break
    return distinct

def max_similar(seq1, seq2):
    similar = 0
    for nt1, nt2 in zip(seq1, seq2):
        if nt1 != nt2:
            break
        similar += 1
    return similar

def create_random_breakpoint(genome, num_inserted, adjacent_length, required_homology=None):
    while True:
        chr1, str1, pos1 = random_chromosome_position(genome, adjacent_length)
        chr2, str2, pos2 = random_chromosome_position(genome, adjacent_length)
        downstream1 = retrieve_downstream_sequence(genome, chr1, str1, pos1, adjacent_length)
        downstream2 = retrieve_downstream_sequence(genome, chr2, str2, pos2, adjacent_length)
        sequence1 = retrieve_upstream_sequence(genome, chr1, str1, pos1, adjacent_length)
        sequence2 = retrieve_upstream_sequence(genome, chr2, str2, pos2, adjacent_length)
        if str1 != '+':
            sequence1 = reverse_complement(sequence1)
            downstream1 = reverse_complement(downstream1)
        if str2 != '-':
            sequence2 = reverse_complement(sequence2)
            downstream2 = reverse_complement(downstream2)
        inserted = create_distinct(downstream1, downstream2, num_inserted)
        homology = max_similar(sequence1[::-1], downstream2[::-1]) + max_similar(sequence2, downstream1)
        if required_homology is not None and homology != required_homology:
            continue
        sequence = sequence1 + inserted + sequence2
        if 'N' in sequence or 'n' in sequence:
            continue
        return chr1, str1, pos1, chr2, str2, pos2, inserted, sequence, homology

def simulate(cfg, sim_info, read_count, sequences_fasta, reads1, reads2, random_reads=True):
    temp_prefix = reads1 + '.dwgsim.temp.prefix'
    random_reads_flag = ''
    if not random_reads:
        random_reads_flag = '-y 0'
    dwgsim_command = cfg.dwgsim_bin + ' -H {7} -z {0} -d {1} -s {2} -N {3} -1 {4} -2 {4} {5} {6}'.format(sim_info['dwgsim_seed'], sim_info['fragment_mean'], sim_info['fragment_stddev'], read_count, sim_info['read_length'], sequences_fasta, temp_prefix, random_reads_flag)
    print dwgsim_command
    dwgsim_retcode = subprocess.call(dwgsim_command.split())
    assert(dwgsim_retcode == 0)
    os.rename(temp_prefix + '.bwa.read1.fastq', reads1)
    os.rename(temp_prefix + '.bwa.read2.fastq', reads2)
    try:
        os.remove(temp_prefix + '.bfast.fastq')
    except OSError:
        pass
    try:
        os.remove(temp_prefix + '.mutations.txt')
    except OSError:
        pass
    try:
        os.remove(temp_prefix + '.mutations.vcf')
    except OSError:
        pass

def create_breakpoints(cfg, sim_info, breakpoints_fasta, breakpoints_info):
    random.seed(int(sim_info['breakpoints_seed']))
    genome = dict(read_sequences(open(cfg.genome_fasta, 'r')))
    with open(breakpoints_fasta, 'w') as fasta, open(breakpoints_info, 'w') as info:
        for idx in range(int(sim_info['num_breakpoints'])):
            chr1, str1, pos1, chr2, str2, pos2, inserted, sequence, homology = create_random_breakpoint(genome, int(sim_info['num_inserted']), int(sim_info['adjacent_length']), int(sim_info['homology']))
            fasta.write('>{0}\n{1}\n'.format(idx, sequence))
            info.write('\t'.join([str(a) for a in (idx, chr1, str1, pos1, chr2, str2, pos2, inserted, homology)]) + '\n')

def create(cfg, sim_info, breakpoints_fasta, breakpoints_info, concordant1, concordant2, discordant1, discordant2):
    create_breakpoints(cfg, sim_info, breakpoints_fasta, breakpoints_info)
    sequences_size = sum([len(seq) for id, seq in read_sequences(open(breakpoints_fasta, 'r'))])
    read_count = int(float(sim_info['coverage']) * sequences_size / float(sim_info['fragment_mean']))
    simulate(cfg, sim_info, read_count, breakpoints_fasta, discordant1, discordant2, False)
    simulate(cfg, sim_info, sim_info['num_concordant'], cfg.genome_fasta, concordant1, concordant2)

