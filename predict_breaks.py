
import collections
import pandas as pd
import numpy as np

import utils.misc


breakpoint_fields = ['cluster_id', 'prediction_id',
                     'chromosome_1', 'strand_1', 'position_1',
                     'chromosome_2', 'strand_2', 'position_2',
                     'inserted']


realignment_fields = ['cluster_id', 'prediction_id',
                      'library_id', 'read_id', 'read_end', 'align_id',
                      'aligned_length', 'template_length', 'score']


score_stats_fields = ['aligned_length', 'expon_lda']


def predict_breaks(clusters_filename, spanning_filename, split_filename, breakpoints_filename):

    # Read all clusters
    fields = ['cluster_id', 'cluster_end', 'lib_id', 'read_id', 'read_end', 'align_id']
    clusters = pd.read_csv(clusters_filename, sep='\t', names=fields)

    reads = clusters[['lib_id', 'read_id']].drop_duplicates()

    def read_filter(df):
        return df.merge(reads, on=['lib_id', 'read_id'], how='inner')

    # Read only spanning reads relevant clusters
    fields = ['lib_id', 'read_id', 'read_end', 'align_id', 'chromosome', 'strand', 'start', 'end', 'score']
    csv_iter = pd.read_csv(spanning_filename, sep='\t', iterator=True, chunksize=1000,
                                              names=fields, converters={'chromosome':str})
    spanning = pd.concat([read_filter(chunk) for chunk in csv_iter])

    span_index_cols = ['lib_id', 'read_id', 'read_end', 'align_id']

    spanning.set_index(span_index_cols, inplace=True)

    # Read only split reads relevant clusters
    fields = ['lib_id', 'read_id', 'read_end',
              'align_id_1', 'chromosome_1', 'strand_1', 'position_1',
              'align_id_2', 'chromosome_2', 'strand_2', 'position_2',
              'inserted', 'score']
    csv_iter = pd.read_csv(split_filename, sep='\t', iterator=True, chunksize=1000,
                                           names=fields, converters={'chromosome':str},
                                           na_values=['.'])
    split = pd.concat([read_filter(chunk) for chunk in csv_iter])

    split_index_cols = ['lib_id', 'read_id', 'align_id_1', 'align_id_2']

    split.set_index(split_index_cols, inplace=True)

    split['inserted'] = split['inserted'].fillna('')

    def flip_split_positions(row):
        if row['flip']:
            row['chromosome_1'], row['chromosome_2'] = row['chromosome_2'], row['chromosome_1']

    predictions = list()

    for cluster_id, cluster_rows in clusters.groupby('cluster_id'):
        
        prediction_id = 0
        
        # Create a table of spanning alignments for this cluster
        cluster_spanning = spanning.merge(cluster_rows, left_index=True, right_on=span_index_cols)
        
        # Predict based on spanning reads
        span_predict_agg = {'chromosome':max, 'strand':max, 'start':min, 'end':max}
        pred = cluster_spanning.groupby('cluster_end').agg(span_predict_agg).reset_index()
        pred['position'] = np.where(pred['strand'] == '+', pred['end'], pred['start'])
        pred = pred.drop(['start', 'end'], axis=1)
        
        # Reformat table
        pred['cluster_id'] = cluster_id
        pred['prediction_id'] = prediction_id
        pred = pred.set_index(['cluster_id', 'prediction_id', 'cluster_end']).unstack()
        pred.columns = [a+'_'+str(b+1) for a, b in pred.columns.values]
        pred.reset_index(inplace=True)
        pred['inserted'] = ''
        
        # Add spanning read prediction
        predictions.append(pred)
        prediction_id += 1
        
        paired = cluster_rows.set_index(['lib_id', 'read_id', 'read_end'])[['cluster_end', 'align_id']].unstack()
        paired.columns = ['cluster_end_1', 'cluster_end_2', 'align_id_1', 'align_id_2']
        
        paired['flip'] = paired['cluster_end_1'] != 0
        paired = paired.drop(['cluster_end_1', 'cluster_end_2'], axis=1)
        paired = paired.reset_index()
        
        cluster_split = split.merge(paired, left_index=True, right_on=split_index_cols)
        
        if len(cluster_split.index) == 0:
            continue
            
        cluster_split.loc[cluster_split['flip'], 'chromosome_1'], cluster_split.loc[cluster_split['flip'], 'chromosome_2'] = \
            cluster_split.loc[cluster_split['flip'], 'chromosome_2'], cluster_split.loc[cluster_split['flip'], 'chromosome_1']
        
        cluster_split.loc[cluster_split['flip'], 'strand_1'], cluster_split.loc[cluster_split['flip'], 'strand_2'] = \
            cluster_split.loc[cluster_split['flip'], 'strand_2'], cluster_split.loc[cluster_split['flip'], 'strand_1']
        
        cluster_split.loc[cluster_split['flip'], 'position_1'], cluster_split.loc[cluster_split['flip'], 'position_2'] = \
            cluster_split.loc[cluster_split['flip'], 'position_2'], cluster_split.loc[cluster_split['flip'], 'position_1']
        
        cluster_split['seed_end'] = np.where(cluster_split['flip'], 1-cluster_split['read_end'], cluster_split['read_end'])
        
        def revcomp_inserted(row):
            if row['seed_end'] == 0:
                return row['inserted']
            else:
                return utils.misc.reverse_complement(row['inserted'])
        
        cluster_split['inserted'] = cluster_split.apply(revcomp_inserted, axis=1)
        
        cluster_split['inslen'] = cluster_split['inserted'].apply(len)
        
        cluster_split.set_index(['position_1', 'position_2', 'inslen'], inplace=True)
        cluster_split = cluster_split.sort_index()
        
        # Calculate highest scoring split
        split_score_sums = cluster_split.groupby(level=[0, 1, 2])['score'].sum()
        split_score_sums.sort(ascending=False)
        
        # Select split alignments for highest scoring split
        pred = cluster_split.loc[split_score_sums.index[0]:split_score_sums.index[0]].reset_index()
        
        # Consensus for inserted sequence
        inserted = np.array([np.array(list(a)) for a in pred['inserted'].values])
        consensus = list()
        for nt_list in inserted.T:
            consensus.append(collections.Counter(nt_list).most_common(1)[0][0])
        consensus = ''.join(consensus)

        # Reformat table
        pred = pred.iloc[0:1].copy()
        pred['cluster_id'] = cluster_id
        pred['prediction_id'] = prediction_id
        pred['inserted'] = consensus
        pred = pred[breakpoint_fields]
        
        predictions.append(pred)
        prediction_id += 1

    if len(predictions) == 0:
        with open(breakpoints_filename, 'w'):
            pass
        return
        
    predictions = pd.concat(predictions, ignore_index=True)

    predictions.loc[predictions['inserted'] == '', 'inserted'] = '.'

    predictions = predictions[breakpoint_fields]
    predictions.to_csv(breakpoints_filename, sep='\t', index=False, header=False)


def calculate_cluster_weights(breakpoints_filename, weights_filename):
    
    epsilon = 0.0001
    itx_distance = 1000000000
    
    breakpoints = pd.read_csv(breakpoints_filename, sep='\t', names=breakpoint_fields)
    
    breakpoints = breakpoints[breakpoints['prediction_id'] == 0]

    breakpoints['distance'] = np.absolute(breakpoints['position_1'] - breakpoints['position_2'])
    breakpoints.loc[breakpoints['chromosome_1'] != breakpoints['chromosome_2'], 'distance'] = itx_distance
    breakpoints['weight'] = 1.0 + epsilon * np.log(breakpoints['distance'])

    breakpoints = breakpoints.sort('cluster_id')
    breakpoints[['cluster_id', 'weight']].to_csv(weights_filename, sep='\t', index=False, header=False)


def calculate_realignment_likelihoods(breakpoints_filename, realignments_filename, score_stats_filename,
                                      likelihoods_filename, match_score, fragment_mean, fragment_stddev):

    match_score = float(match_score)
    fragment_mean = float(fragment_mean)
    fragment_stddev = float(fragment_stddev)

    score_stats = pd.read_csv(score_stats_filename, sep='\t', names=score_stats_fields)

    realignments = pd.read_csv(realignments_filename, sep='\t', names=realignment_fields)

    breakpoints = pd.read_csv(breakpoints_filename, sep='\t', names=breakpoint_fields,
                              converters={'chromosome_1':str, 'chromosome_2':str},
                              na_values=['.'])

    breakpoints['inserted'] = breakpoints['inserted'].fillna('')

    breakpoints['inslen'] = breakpoints['inserted'].apply(len)

    data = realignments.merge(breakpoints[['cluster_id', 'prediction_id', 'inslen']],
                              on=['cluster_id', 'prediction_id'])

    data = data.drop_duplicates(['cluster_id', 'prediction_id', 'read_id', 'read_end'])

    data = data.merge(score_stats, on='aligned_length')

    data['max_score'] = match_score * data['aligned_length']
    data['score_diff'] = data['max_score'] - data['score']
    data['score_log_likelihood'] = np.log(data['expon_lda']) - \
                                           data['expon_lda'] * data['score_diff']

    agg_f = {'score_log_likelihood':sum,
             'template_length':sum,
             'inslen':max}
    data = data.groupby(['cluster_id', 'prediction_id', 'read_id']).agg(agg_f)
    data['template_length'] += data['inslen']

    constant = 1. / ((2 * np.pi)**0.5 * fragment_stddev)
    data['length_log_likelihood'] = -np.log(constant) - \
                                            np.square(data['template_length'] - fragment_mean) / (2. * fragment_stddev**2)

    data['log_likelihood'] = data['score_log_likelihood'] + data['length_log_likelihood']

    data = data['log_likelihood'].reset_index()

    data.to_csv(likelihoods_filename, sep='\t', index=False)



