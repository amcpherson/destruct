import string
import numpy as np


def reverse_complement(sequence):
    return sequence[::-1].translate(str.maketrans('ACTGactg','TGACtgac'))

def column_flip(df, cond, col_1, col_2):
    df.loc[cond, col_1], df.loc[cond, col_2] = df.loc[cond, col_2], df.loc[cond, col_1]


def calculate_offset(strand, offset):
    if (strand == "+"):
        direction = 1
    else:
        direction = -1
    return offset * direction


def calculate_forward_homology(chromosome, strand, position, genome, maxOffset, flip=False):
    if flip:
        idx1 = 1
    else:
        idx1 = 0
    idx2 = 1 - idx1

    chr1 = genome[chromosome[idx1]]
    chr2 = genome[chromosome[idx2]]

    homology = 0
    for offset in range(1, maxOffset + 1):
        nt1 = chr1[position[idx1] + calculate_offset(strand[idx1], offset) - 1]
        nt2 = chr2[position[idx2] + calculate_offset(strand[idx2], 1 - offset) - 1]

        if (strand[idx1] != '+'):
            nt1 = reverse_complement(nt1)
        
        if (strand[idx2] != '-'):
            nt2 = reverse_complement(nt2)

        if (nt1 != nt2):
            break
        
        homology = offset
    
    return homology


def homology_consistent_breakpoint(chromosome, strand, position, genome, maxOffset):
    maxOffsetA = calculate_forward_homology(chromosome, strand, position, genome, maxOffset, False)
    maxOffsetB = calculate_forward_homology(chromosome, strand, position, genome, maxOffset, True)

    # Ensure that the same breakpoint is selected among the multiple
    # breakpoints possible when there is breakpoint homology.  Always
    # select the breakpoint for which the minimum of the two breakend
    # positions is minimal

    positionA1 = position[0] + calculate_offset(strand[0], maxOffsetA)
    positionA2 = position[1] + calculate_offset(strand[1], -maxOffsetA)

    positionB1 = position[0] + calculate_offset(strand[0], -maxOffsetB)
    positionB2 = position[1] + calculate_offset(strand[1], maxOffsetB)

    if (min(positionA1, positionA2) < min(positionB1, positionB2)):
        position[0] = positionA1
        position[1] = positionA2
    else:
        position[0] = positionB1
        position[1] = positionB2

    homology = maxOffsetA + maxOffsetB

    return homology


def normalize_breakpoint(chr1, str1, pos1, chr2, str2, pos2, genome, max_offset=100):
    chromosome = [chr1, chr2]
    strand = [str1, str2]
    position = [pos1, pos2]

    homology = homology_consistent_breakpoint(
        chromosome, strand, position, genome, max_offset)

    return position[0], position[1], homology

