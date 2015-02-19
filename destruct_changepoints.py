
import argparse    
import pandas as pd

argparser = argparse.ArgumentParser()

argparser.add_argument('input_breakpoint', help='destruct breakpoint file')

argparser.add_argument('input_breakpoint_library', help='destruct breakpoint library file')

argparser.add_argument('control_id', help='identifier of control sample')

argparser.add_argument('output_changepoints', help='prefix for plot files')

args = argparser.parse_args()

brk = pd.read_csv(args.input_breakpoint, sep='\t',
                  converters={'chromosome_1':str, 'chromosome_2':str})

brklib = pd.read_csv(args.input_breakpoint_library, sep='\t')


# Add is_normal column
brklib['is_normal'] = brklib['library'] == args.control_id


# Mark as germline any prediction with nonzero normal reads

is_germline = brklib[brklib['num_reads'] > 0].groupby('prediction_id')['is_normal'].any()

brk.set_index('prediction_id', inplace=True)
brk['is_germline'] = is_germline
brk.reset_index(inplace=True)


# Mark as a filtered any breakpoint that is common or germline

brk['is_filtered'] = brk['is_germline']


# Get a list of breakends and their filtered status

def get_brkend(brk, side, data_cols):
    cols = ['chromosome', 'strand', 'position']
    side_cols = [a+'_'+side for a in cols]
    brkend = brk[['prediction_id']+side_cols+data_cols]
    brkend = brkend.rename(columns=dict(zip(side_cols, cols)))
    brkend['side'] = side
    return brkend

brkend = pd.concat([get_brkend(brk, '1', ['is_filtered']),
                    get_brkend(brk, '2', ['is_filtered'])], ignore_index=True)

# Calculate the minimum distance to the nearest filtered breakend for all breakpoints

def calculate_dist_filtered(data):

    data = data.sort('position')

    data['nearest_left'] = None
    data['nearest_right'] = None

    data.loc[data['is_filtered'], 'nearest_left'] = data.loc[data['is_filtered'], 'position']
    data.loc[data['is_filtered'], 'nearest_right'] = data.loc[data['is_filtered'], 'position']

    data['nearest_left'] = data['nearest_left'].fillna(method='ffill')
    data['nearest_right'] = data['nearest_right'].fillna(method='bfill')

    data['nearest_left'] = data['position'] - data['nearest_left']
    data['nearest_right'] = data['nearest_right'] - data['position']

    data['dist_filtered'] = data[['nearest_left', 'nearest_right']].min(axis=1)
    
    return data['dist_filtered']

dist_filtered = brkend.set_index(['prediction_id', 'side'])\
                      .groupby(['chromosome', 'strand'])\
                      .apply(calculate_dist_filtered)\
                      .reset_index(level=[0,1], drop=True)\
                      .unstack(level=1)\
                      .min(axis=1)

brk.set_index('prediction_id', inplace=True)
brk['dist_filtered'] = dist_filtered
brk.reset_index(inplace=True)

chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
               '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22',
               'X', 'Y']

# Filter breakpoints

brk = brk.loc[(~brk['is_filtered']) &
              (brk['dist_filtered'] > 50) &
              (brk['num_reads'] >= 2) &
              (brk['num_split'] >= 2) &
              (brk['log_likelihood'] > -20.) &
              (brk['template_length_min'] > 120) &
              (brk['chromosome_1'].isin(chromosomes)) &
              (brk['chromosome_2'].isin(chromosomes))]


# Filter brklib table

brklib = brklib.merge(brk[['prediction_id']].drop_duplicates())


# Output changepoints
changepoints = pd.concat([brk[['chromosome_1', 'position_1']].rename(columns=lambda a: a[:-2]), 
                          brk[['chromosome_2', 'position_2']].rename(columns=lambda a: a[:-2])], ignore_index=True)
    
changepoints['position'] = changepoints['position'].astype(int)

changepoints.to_csv(args.output_changepoints, sep='\t', header=False, index=False)

