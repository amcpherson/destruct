import pandas as pd
import numpy as np


breakpoint_fields = ['cluster_id', 'prediction_id',
                     'chromosome_1', 'strand_1', 'position_1',
                     'chromosome_2', 'strand_2', 'position_2']

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

    # Read only spanning reads relevant clusters
    fields = ['lib_id', 'read_id', 'read_end',
              'align_id_1', 'chromosome_1', 'strand_1', 'position_1',
              'align_id_2', 'chromosome_2', 'strand_2', 'position_2',
              'inserted', 'score']
    csv_iter = pd.read_csv(split_filename, sep='\t', iterator=True, chunksize=1000,
                                           names=fields, converters={'chromosome':str})
    split = pd.concat([read_filter(chunk) for chunk in csv_iter])

    split_index_cols = ['lib_id', 'read_id', 'align_id_1', 'align_id_2']

    split.set_index(split_index_cols, inplace=True)

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
        
        predictions.append(pred)
        prediction_id += 1
        
        paired = cluster_rows.set_index(['lib_id', 'read_id', 'cluster_end'])[['read_end', 'align_id']].unstack()
        paired.columns = ['read_end_1', 'read_end_2', 'align_id_1', 'align_id_2']
        
        paired['flip'] = paired['read_end_1'] != 0
        paired = paired.drop(['read_end_1', 'read_end_2'], axis=1)
        paired = paired.reset_index()
        
        cluster_split = split.merge(paired, left_index=True, right_on=split_index_cols)
        
        cluster_split.loc[cluster_split['flip'], 'chromosome_1'], cluster_split.loc[cluster_split['flip'], 'chromosome_2'] =         cluster_split.loc[cluster_split['flip'], 'chromosome_2'], cluster_split.loc[cluster_split['flip'], 'chromosome_1']
        
        cluster_split.loc[cluster_split['flip'], 'strand_1'], cluster_split.loc[cluster_split['flip'], 'strand_2'] =         cluster_split.loc[cluster_split['flip'], 'strand_2'], cluster_split.loc[cluster_split['flip'], 'strand_1']
        
        cluster_split.loc[cluster_split['flip'], 'position_1'], cluster_split.loc[cluster_split['flip'], 'position_2'] =         cluster_split.loc[cluster_split['flip'], 'position_2'], cluster_split.loc[cluster_split['flip'], 'position_1']
        
        if len(cluster_split.index) == 0:
            continue
            
        cluster_split = cluster_split.drop('flip', axis=1)
        
        cluster_split.set_index(['position_1', 'position_2'], inplace=True)
        cluster_split['score_sum'] = cluster_split.groupby(level=[0, 1])['score'].sum()
        cluster_split.reset_index(inplace=True)
        
        pred = cluster_split[cluster_split['score_sum'] == cluster_split['score_sum'].max()].iloc[0:1].copy()
        pred['cluster_id'] = cluster_id
        pred['prediction_id'] = prediction_id
        pred = pred[breakpoint_fields]
        
        predictions.append(pred)
        prediction_id += 1

    if len(predictions) == 0:
        with open(breakpoints_filename, 'w'):
            pass
        return
        
    predictions = pd.concat(predictions, ignore_index=True)

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


