import random
import os
import subprocess
import pandas as pd

import destruct.utils.seq
import destruct.utils.misc


def random_chromosome_position(genome, dist_to_end):
    """ Generate random chromosome, strand position.

     > Note that positions are 1-based
    """

    genome_length = sum([len(chr) for chr in genome.values()])
    while True:
        position = random.randint(1, genome_length + 1)
        for id, chr in genome.items():
            if position <= len(chr):
                if position <= dist_to_end or len(chr) - position < dist_to_end:
                    break
                return id, random.choice(['+', '-']), position
            position -= len(chr)


def retrieve_upstream_sequence(genome, chromosome, strand, position, length):
    """ Retrieve sequence upstream of a position, including that position
    """

    if strand == '+':
        end = position
        start = position - length + 1
    else:
        start = position
        end = position + length - 1

    # Positions assumed 1-based, closed region
    start -= 1

    return genome[chromosome][start:end]


def retrieve_downstream_sequence(genome, chromosome, strand, position, length):
    """ Retrieve sequence downstream of a position, not including that position
    """

    if strand == '+':
        start = position + 1
        end = start + length - 1
    else:
        end = position - 1
        start = end - length + 1

    # Positions assumed 1-based, closed region
    start -= 1

    return genome[chromosome][start:end]


def create_distinct(sequence1, sequence2, length):
    """ Create an insertion sequence distinct from reference genome sequence
    """

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
    """ Calculate the maximum number of identical nucleotides
    """

    similar = 0
    for nt1, nt2 in zip(seq1, seq2):
        if nt1 != nt2:
            break
        similar += 1
    return similar


def create_random_breakpoint(genome, adjacent_length, num_inserted, required_homology=None):
    """ Create a random breakpoint with the required parameters
    """

    while True:
        chr1, str1, pos1 = random_chromosome_position(genome, adjacent_length)
        chr2, str2, pos2 = random_chromosome_position(genome, adjacent_length)
        downstream1 = retrieve_downstream_sequence(genome, chr1, str1, pos1, adjacent_length)
        downstream2 = retrieve_downstream_sequence(genome, chr2, str2, pos2, adjacent_length)
        sequence1 = retrieve_upstream_sequence(genome, chr1, str1, pos1, adjacent_length)
        sequence2 = retrieve_upstream_sequence(genome, chr2, str2, pos2, adjacent_length)
        if str1 != '+':
            sequence1 = destruct.utils.misc.reverse_complement(sequence1)
            downstream1 = destruct.utils.misc.reverse_complement(downstream1)
        if str2 != '-':
            sequence2 = destruct.utils.misc.reverse_complement(sequence2)
            downstream2 = destruct.utils.misc.reverse_complement(downstream2)
        inserted = create_distinct(downstream1, downstream2, num_inserted)
        homology = max_similar(sequence1[::-1], downstream2[::-1]) + max_similar(sequence2, downstream1)
        if required_homology is not None and homology != required_homology:
            continue
        sequence = sequence1 + inserted + sequence2
        if 'N' in sequence or 'n' in sequence:
            continue
        if inserted == '':
            norm_pos1, norm_pos2, homology_check = destruct.utils.misc.normalize_breakpoint(chr1, str1, pos1, chr2, str2, pos2, genome)
            if homology_check != homology:
                raise Exception(
                    '''homology {} doesnt match expected homology {} for 
                    breakpoint {} {} {}, {} {} {}'''.format(
                        homology_check, homology,
                        chr1, str1, pos1, chr2, str2, pos2))
            pos1, pos2 = norm_pos1, norm_pos2
        return chr1, str1, pos1, chr2, str2, pos2, inserted, sequence, homology

def simulate(sim_info, read_count, sequences_fasta, reads1, reads2, random_reads=True):
    """ Simulate reads
    """

    temp_prefix = reads1 + '.dwgsim.temp.prefix'
    dwgsim_command = ['dwgsim']
    dwgsim_command += ['-z', sim_info['dwgsim_seed']]
    dwgsim_command += ['-d', sim_info['fragment_mean']]
    dwgsim_command += ['-s', sim_info['fragment_stddev']]
    dwgsim_command += ['-N', read_count]
    dwgsim_command += ['-1', sim_info['read_length']]
    dwgsim_command += ['-2', sim_info['read_length']]
    dwgsim_command += ['-H']
    if not random_reads:
        dwgsim_command += ['-y', '0']
    if sim_info.get('perfect_reads', '') == 'yes':
        dwgsim_command += ['-e', '0']
        dwgsim_command += ['-E', '0']
        dwgsim_command += ['-r', '0']
        dwgsim_command += ['-R', '0']
    dwgsim_command += [sequences_fasta]
    dwgsim_command += [temp_prefix]
    dwgsim_command = [str(a) for a in dwgsim_command]
    print (' '.join(dwgsim_command))
    dwgsim_retcode = subprocess.call(dwgsim_command)
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

def create_breakpoints(sim_info, genome_fasta, breakpoints_fasta, breakpoints_info):
    """ Create a set of simulated breakpoints
    """

    random.seed(int(sim_info['breakpoints_seed']))
    genome = dict(destruct.utils.seq.read_sequences(open(genome_fasta, 'r')))
    info_table = []
    with open(breakpoints_fasta, 'w') as fasta:
        for idx in range(int(sim_info['num_breakpoints'])):
            adjacent_length = int(sim_info['adjacent_length'])
            if sim_info['random_break_features']:
                homology = random.randint(0, int(sim_info['homology']))
                if homology == 0 and random.random() < 0.1:
                    num_inserted = random.randint(0, int(sim_info['num_inserted']))
                else:
                    num_inserted = 0
            else:
                num_inserted = int(sim_info['num_inserted'])
                homology = int(sim_info['homology'])
            chr1, str1, pos1, chr2, str2, pos2, inserted, sequence, homology = create_random_breakpoint(
                genome, adjacent_length, num_inserted, homology)
            fasta.write('>{0}\n{1}\n'.format(idx, sequence))
            info_table.append(
                {
                    'breakpoint_id': idx,
                    'chromosome_1': chr1,
                    'strand_1': str1,
                    'position_1': pos1,
                    'chromosome_2': chr2,
                    'strand_2': str2,
                    'position_2': pos2,
                    'inserted': inserted,
                    'homology': homology,
                }
            )
    info_table = pd.DataFrame(info_table)
    info_table.to_csv(breakpoints_info, sep='\t', index=False)

def create(sim_info, genome_fasta, breakpoints_fasta, breakpoints_info, concordant1, concordant2, discordant1, discordant2):
    """ Create simulated breakpoints and simulate reads from those breakpoints
    """

    create_breakpoints(sim_info, genome_fasta, breakpoints_fasta, breakpoints_info)
    sequences_size = sum([len(seq) for id, seq in destruct.utils.seq.read_sequences(open(breakpoints_fasta, 'r'))])
    read_count = int(float(sim_info['coverage']) * sequences_size / float(sim_info['fragment_mean']))
    simulate(sim_info, read_count, breakpoints_fasta, discordant1, discordant2, False)
    simulate(sim_info, sim_info['num_concordant'], genome_fasta, concordant1, concordant2)
