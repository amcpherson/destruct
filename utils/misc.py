import string
import numpy as np

def reverse_complement(sequence):
    return sequence[::-1].translate(string.maketrans('ACTGactg','TGACtgac'))

def column_flip(df, cond, col_1, col_2):
    df.loc[cond, col_1], df.loc[cond, col_2] = df.loc[cond, col_2], df.loc[cond, col_1]

def normalize_breakpoints(data):

    # Calculate advancement due to homology based on strand
    data['advance_1'] = np.where(data['strand_1'] == "+", data['homology'], -data['homology'])
    data['advance_2'] = np.where(data['strand_2'] == "+", data['homology'], -data['homology'])

    # Add to the first and subtract from the second to create alternate breakends
    data['position_alt_1'] = data['position_1'] + data['advance_1']
    data['position_alt_2'] = data['position_2'] + data['advance_2']

    # Flip based on lexical ordering of chromosome, strand, position
    # Ensure that the positions are such that the alternate breakpoint
    # can be obtained by adding the homology to position_1 and subtracting
    # the homology from position_2
    data['flip'] = False
    data['flip'] |= (data['chromosome_2'] < data['chromosome_1'])
    data['flip'] |= ((data['chromosome_2'] == data['chromosome_1']) & \
                     (data['strand_2'] < data['strand_1']))
    data['flip'] |= ((data['chromosome_2'] == data['chromosome_1']) & \
                     (data['strand_2'] == data['strand_1']) & \
                     (data['position_alt_2'] < data['position_1']))

    # Flip columns to make them consistent across clusters
    column_flip(data, data['flip'], 'chromosome_1', 'chromosome_2')
    column_flip(data, data['flip'], 'strand_1', 'strand_2')
    data.loc[data['flip'], 'position_1'] = data.loc[data['flip'], 'position_alt_2']
    data.loc[data['flip'], 'position_2'] = data.loc[data['flip'], 'position_alt_1']

    data.drop(['flip', 'advance_1', 'advance_2', 'position_alt_1', 'position_alt_2'], axis=1, inplace=True)
        
    return data
