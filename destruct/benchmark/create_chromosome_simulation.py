import argparse
import ConfigParser
import random
import string
import os
import subprocess
import sys

import utils.seq
import utils.misc


class SequenceRegion(object):
    def __init__(self, chromosome, strand, start, end):
        self.chromosome = chromosome
        self.strand = strand
        self.start = start
        self.end = end
    def __repr__(self):
        return 'SequenceRegion({0},{1},{2},{3})'.format(repr(self.chromosome), repr(self.strand), repr(self.start), repr(self.end))
    def __eq__(self, other):
        return self.chromosome == other.chromosome and self.strand == other.strand and self.start == other.start and self.end == other.end
    @property
    def length(self):
        return self.end - self.start + 1
    def subregion(self, local_start, local_end):
        return SequenceRegion(self.chromosome, self.strand, self.start + local_start - 1, self.start + local_end - 1)
    def invert(self):
        strand = ('+', '-')[self.strand == '+']
        assert strand != self.strand
        return SequenceRegion(self.chromosome, strand, self.start, self.end)

class TumourChromosome(object):
    def __init__(self, regions, genome):
        self.regions = regions
        self.length = sum([a.end - a.start + 1 for a in regions])
        self.genome = genome
    def local_start_ends(self):
        """
        generate a list of regions, in addition to the start and end
        position of the region in the tumour chromosome
        """
        offset = 0
        for region in self.regions:
            yield region, offset + 1, offset + region.length
            offset += region.length
    def reference_position(self, position):
        """ given a position in the tumour chromosome, identify the position
        in the reference genome
        """
        assert position <= self.length
        offset = 0
        for region in self.regions:
            if position <= offset + region.length:
                region_position = position - offset
                if region.strand == '-':
                    region_position = region.length - region_position + 1
                return region.chromosome, region.strand, region.start + region_position - 1
            offset += region.length
    def generate_spaced_intervals(self, num, size_func, spacing=1000):
        """ generate a set of non overlapping intervals in the tumour
        chromosome such that there is at least spacing nucleotides between
        intervals and at least spacing nucleotides between interval boundaries
        and existing distinct sequence regions that make up the tumour
        chromosome
        """
        intervals = list()
        while len(intervals) < num:
            size = size_func()
            start = random.randint(1 + spacing, self.length - size - spacing - 1)
            end = start + size - 1
            overlapping = False
            for info in intervals:
                if start <= info[1] + spacing and end >= info[0] - spacing:
                    overlapping = True
                    break
            for region, local_start, local_end in self.local_start_ends():
                if abs(start - local_start) < spacing or abs(start - local_end) < spacing or abs(end - local_start) < spacing or abs(end - local_end) < spacing:
                    overlapping = True
                    break
            if overlapping:
                continue
            unmappable = False
            for chromosome_position in (start, end):
                chromosome, strand, position = self.reference_position(chromosome_position)
                if self.genome[chromosome].count('N', position - 100, position + 100) > 0:
                    unmappable = True
            if unmappable:
                continue
            intervals.append((start, end))
        return sorted(intervals)
    def split_regions(self, intervals):
        """ split the current tumour chromosome regions at the boundaries
        of the given set of intervals
        """
        for interval in intervals:
            new_regions = list()
            for region, local_start, local_end in self.local_start_ends():
                start_ends = list([1, region.length])
                if interval[0] > local_start and interval[0] < local_end:
                    start_ends.append(interval[0] - local_start)
                    start_ends.append(interval[0] - local_start + 1)
                if interval[1] > local_start and interval[1] < local_end:
                    start_ends.append(interval[1] - local_start + 1)
                    start_ends.append(interval[1] - local_start + 2)
                start_ends = sorted(start_ends)
                for local_start, local_end in zip(start_ends[::2], start_ends[1::2]):
                    new_regions.append(region.subregion(local_start, local_end))
            self.regions = sorted(new_regions, key=lambda a: a.start)
    def merge_regions(self):
        """ merge adjacent tumour chromosome regions if those regions are adjacent
        and consecutive in the reference genome
        """
        new_regions = list()
        idx1 = 0
        while idx1 < len(self.regions):
            idx2 = idx1 + 1
            while idx2 < len(self.regions):
                if self.regions[idx1].chromosome != self.regions[idx2].chromosome:
                    break
                if self.regions[idx1].strand != self.regions[idx2].strand:
                    break
                if self.regions[idx2-1].end != self.regions[idx2].start - 1:
                    break
                idx2 += 1
            idx2 -= 1
            new_regions.append(SequenceRegion(self.regions[idx1].chromosome, self.regions[idx1].strand, self.regions[idx1].start, self.regions[idx2].end))
            idx1 = idx2 + 1
        self.regions = new_regions
        assert self.length == sum([a.end - a.start + 1 for a in new_regions])
    def simulate_deletions(self, num, size_func, spacing=1000):
        """ simulate deletions of size given by size_func spaced by spacing
        """
        deleted = self.generate_spaced_intervals(num, size_func, spacing)
        self.split_regions(deleted)
        new_regions = list()
        for region, local_start, local_end in self.local_start_ends():
            if len(deleted) == 0 or local_end < deleted[0][0]:
                new_regions.append(region)
            elif local_end <= deleted[0][1]:
                pass
            elif local_end > deleted[0][1]:
                deleted = deleted[1:]
                new_regions.append(region)
        self.regions = new_regions
        self.length = sum([a.end - a.start + 1 for a in new_regions])
        self.merge_regions()
    def simulate_duplications(self, num, size_func, spacing=1000):
        """ simulate duplications of size given by size_func spaced by spacing
        """
        duplicated = self.generate_spaced_intervals(num, size_func, spacing)
        self.split_regions(duplicated)
        new_regions = list()
        dup_regions = list()
        for region, local_start, local_end in self.local_start_ends():
            if len(duplicated) == 0 or local_end < duplicated[0][0]:
                new_regions.append(region)
            elif local_end <= duplicated[0][1]:
                dup_regions.append(region)
            elif local_end > duplicated[0][1]:
                new_regions.extend(dup_regions)
                new_regions.extend(dup_regions)
                dup_regions = list()
                duplicated = duplicated[1:]
                new_regions.append(region)
        self.regions = new_regions
        self.length = sum([a.end - a.start + 1 for a in new_regions])
        self.merge_regions()
    def simulate_inversions(self, num, size_func, spacing=1000):
        """ simulate inversions of size given by size_func spaced by spacing
        """
        inverted = self.generate_spaced_intervals(num, size_func, spacing)
        self.split_regions(inverted)
        new_regions = list()
        inv_regions = list()
        for region, local_start, local_end in self.local_start_ends():
            if len(inverted) == 0 or local_end < inverted[0][0]:
                new_regions.append(region)
            elif local_end <= inverted[0][1]:
                inv_regions.append(region)
            elif local_end > inverted[0][1]:
                new_regions.extend([a.invert() for a in reversed(inv_regions)])
                inv_regions = list()
                inverted = inverted[1:]
                new_regions.append(region)
        self.regions = new_regions
        self.length = sum([a.end - a.start + 1 for a in new_regions])
        self.merge_regions()
    def create_sequence(self, genome):
        sequence = ''
        for region in self.regions:
            region_sequence = genome[region.chromosome][region.start-1:region.end]
            if region.strand == '-':
                region_sequence = utils.misc.reverse_complement(region_sequence)
            sequence += region_sequence
        return sequence
    @property
    def breakpoints(self):
        for id in range(len(self.regions)-1):
            region1 = self.regions[id]
            region2 = self.regions[id+1]
            if region1.strand == '+':
                pos1 = region1.end
                str1 = '+'
            else:
                pos1 = region1.start
                str1 = '-'
            if region2.strand == '+':
                pos2 = region2.start
                str2 = '-'
            else:
                pos2 = region2.end
                str2 = '+'
            yield id, region1.chromosome, str1, pos1, region2.chromosome, str2, pos2

def random_chimera(fasta, chromosome, num_chimera, fragment_length):
    id = 0
    chimera_length = 0
    while id < num_chimera:
        chimera = ''
        for side in (0, 1):
            start = random.randint(0, len(chromosome) - fragment_length - 1)
            end = start + fragment_length
            sequence = chromosome[start:end]
            if random.choice([True, False]):
                sequence = utils.misc.reverse_complement(sequence)
            chimera += sequence
        if chimera.count('N') > 0:
            continue
        fasta.write('>r{0}\n{1}\n'.format(id, chimera))
        id += 1
        chimera_length += len(chimera)
    return chimera_length

def write_fasta(fasta, id, sequence):
    nt_per_line = 80
    fasta.write('>' + str(id) + '\n')
    for start in range(0, len(sequence), nt_per_line):
        fasta.write(sequence[start:start+nt_per_line] + '\n')

def simulate(cfg, sim_info, read_count, sequences_fasta, reads1, reads2, random_reads=True):
    print ('Simulating {0} reads from {1}'.format(read_count, sequences_fasta))
    temp_prefix = reads1 + '.dwgsim.temp.prefix'
    random_reads_flag = ''
    if not random_reads:
        random_reads_flag = '-y 0'
    dwgsim_command = cfg.dwgsim_bin + ' -H {7} -z {0} -d {1} -s {2} -N {3} -1 {4} -2 {4} {5} {6}'.format(sim_info['dwgsim_seed'], sim_info['fragment_mean'], sim_info['fragment_stddev'], read_count, sim_info['read_length'], sequences_fasta, temp_prefix, random_reads_flag)
    print (dwgsim_command)
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

def create(cfg, sim_info, chromosomes_info, breakpoints_info, reads_fastq_1, reads_fastq_2, temps_prefix):
    random.seed(int(sim_info['rearrangements_seed']))
    genome = dict(utils.seq.read_sequences(open(cfg.genome_fasta, 'r')))
    simulation_chromosome = sim_info['simulation_chromosome']
    # Create simulated chromosomes by applying random operations
    tumour_proportions = [float(a) for a in sim_info['chromosome_proportions'].split()]
    simulated_chromosomes = list()
    for chromosome_id, proportion in enumerate(tumour_proportions):
        simulated_chromosome = TumourChromosome([SequenceRegion(simulation_chromosome, '+', 1, len(genome[simulation_chromosome]))], genome)
        simulated_chromosome.simulate_deletions(int(sim_info['deletion_count']), lambda: random.randint(*[int(a) for a in sim_info['deletion_range'].split()]))
        simulated_chromosome.simulate_duplications(int(sim_info['duplication_count']), lambda: random.randint(*[int(a) for a in sim_info['duplication_range'].split()]))
        simulated_chromosome.simulate_inversions(int(sim_info['inversion_count']), lambda: random.randint(*[int(a) for a in sim_info['inversion_range'].split()]))
        simulated_chromosomes.append(simulated_chromosome)
    # Write chromosome regions info
    with open(chromosomes_info, 'w') as info:
        for chromosome_id, simulated_chromosome in enumerate(simulated_chromosomes):
            for region in simulated_chromosome.regions:
                info.write('\t'.join([str(a) for a in (chromosome_id, region.chromosome, region.strand, region.start, region.end)]) + '\n')
    # Write breakpoint positions info
    with open(breakpoints_info, 'w') as info:
        for chromosome_id, simulated_chromosome in enumerate(simulated_chromosomes):
            for break_info in simulated_chromosome.breakpoints:
                info.write('\t'.join([str(a) for a in (chromosome_id, tumour_proportions[chromosome_id]) + break_info]) + '\n')
    create_sequences(cfg, sim_info, genome, simulated_chromosomes, tumour_proportions, reads_fastq_1, reads_fastq_2, temps_prefix)

def create_sequences(cfg, sim_info, genome, chromosomes, proportions, fastq1_filename, fastq2_filename, temps_prefix):
    # Generate fastas from simulated chromosomes
    chromosome_fastas = dict()
    for chromosome_id in range(len(chromosomes)):
        chromosome_fastas[chromosome_id] = temps_prefix + str(chromosome_id) + '.fa'
        with open(chromosome_fastas[chromosome_id], 'w') as fasta:
            write_fasta(fasta, 't'+str(chromosome_id), chromosomes[chromosome_id].create_sequence(genome))
    # Simulate reads from each tumour chromosome
    fastqs1 = list()
    fastqs2 = list()
    for chromosome_id in range(len(chromosomes)):
        scaled_coverage = proportions[chromosome_id] * float(sim_info['coverage'])
        scaled_read_count = int(float(scaled_coverage) * chromosomes[chromosome_id].length / float(sim_info['fragment_mean']))
        temp_fastq1_filename = temps_prefix + str(chromosome_id) + '.1.fastq'
        temp_fastq2_filename = temps_prefix + str(chromosome_id) + '.2.fastq'
        simulate(cfg, sim_info, scaled_read_count, chromosome_fastas[chromosome_id], temp_fastq1_filename, temp_fastq2_filename)
        fastqs1.append(temp_fastq1_filename)
        fastqs2.append(temp_fastq2_filename)
    # Simulate random chimeric fragments
    false_chimera_fasta = temps_prefix + 'false_chimera.fa'
    with open(false_chimera_fasta, 'w') as fasta:
        chimera_length = random_chimera(fasta, genome[sim_info['simulation_chromosome']], int(sim_info['num_false_chimera']), int(sim_info['fragment_mean']))
    false_chimera_fastq1 = temps_prefix + 'f.1.fastq'
    false_chimera_fastq2 = temps_prefix + 'f.2.fastq'
    fastqs1.append(false_chimera_fastq1)
    fastqs2.append(false_chimera_fastq2)
    false_chimera_read_count = int(float(sim_info['false_chimera_coverage']) * chimera_length / float(sim_info['fragment_mean']))
    simulate(cfg, sim_info, false_chimera_read_count, false_chimera_fasta, false_chimera_fastq1, false_chimera_fastq2, False)
    # Concatenate fastqs
    with open(fastq1_filename, 'w') as fastq1:
        subprocess.check_call(['cat'] + fastqs1, stdout=fastq1)
    with open(fastq2_filename, 'w') as fastq2:
        subprocess.check_call(['cat'] + fastqs2, stdout=fastq2)
    for fastq in fastqs1 + fastqs2:
        os.remove(fastq)

