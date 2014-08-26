
import collections
import pandas as pd
import numpy as np
import scipy
import scipy.stats

import utils.misc


cluster_fields = ['cluster_id', 'cluster_end',
                  'library_id', 'read_id', 'read_end', 'align_id']


spanning_fields = ['library_id', 'read_id', 'read_end', 'align_id',
                   'chromosome', 'strand', 'start', 'end', 'score']


split_fields = ['library_id', 'read_id', 'read_end',
                'align_id_1', 'chromosome_1', 'strand_1', 'position_1',
                'align_id_2', 'chromosome_2', 'strand_2', 'position_2',
                'homology', 'inserted', 'score']


breakpoint_fields = ['cluster_id', 'prediction_id',
                     'chromosome_1', 'strand_1', 'position_1',
                     'chromosome_2', 'strand_2', 'position_2',
                     'homology', 'count', 'inserted']


realignment_fields = ['cluster_id', 'prediction_id', 'cluster_end',
                      'library_id', 'read_id', 'read_end', 'align_id',
                      'aligned_length', 'template_length', 'score']


score_stats_fields = ['aligned_length', 'expon_lda']


likelihoods_fields = ['cluster_id', 'prediction_id',
                      'library_id', 'read_id',
                      'read_end_1', 'read_end_2',
                      'aligned_length_1', 'aligned_length_2',
                      'template_length_1', 'template_length_2',
                      'score_1', 'score_2',
                      'score_log_likelihood_1', 'score_log_likelihood_2',
                      'score_log_cdf_1', 'score_log_cdf_2',
                      'inslen', 'template_length', 'length_z_score',
                      'length_log_likelihood', 'length_log_cdf',
                      'log_likelihood', 'log_cdf']


def predict_breaks_spanning(clusters, spanning):

    # Create a table of spanning alignments for this cluster
    data = spanning.merge(clusters, on=['library_id', 'read_id', 'read_end', 'align_id'])

    # Predict based on spanning reads
    agg_f = {'chromosome':max, 'strand':max, 'start':min, 'end':max}
    data = data.groupby(['cluster_id', 'cluster_end']).agg(agg_f).reset_index()
    data['position'] = np.where(data['strand'] == '+', data['end'], data['start'])
    data = data.drop(['start', 'end'], axis=1)

    # Reformat table
    data = data.set_index(['cluster_id', 'cluster_end']).unstack()
    data.columns = [a+'_'+str(b+1) for a, b in data.columns.values]
    data.reset_index(inplace=True)
    data['count'] = 0
    data['homology'] = 0
    data['inserted'] = ''

    return data


def predict_breaks_split(clusters, split):

    # Unstack clusters file by read end, in preparation for
    # merge with split alignments
    paired = clusters.set_index(['cluster_id', 'library_id', 'read_id', 'read_end'])[['cluster_end', 'align_id']].unstack()
    paired.columns = ['cluster_end_1', 'cluster_end_2', 'align_id_1', 'align_id_2']

    # Keep track of which reads have ends flipped relative to clusters
    paired['flip'] = paired['cluster_end_1'] != 0
    paired = paired.drop(['cluster_end_1', 'cluster_end_2'], axis=1)
    paired = paired.reset_index()

    # Merge split read information into clusters
    data = split.merge(paired, on=['library_id', 'read_id', 'align_id_1', 'align_id_2'])
    data['position_1'] = data['position_1'].astype(int)
    data['position_2'] = data['position_2'].astype(int)

    # Flip columns to make them consistent across clusters
    utils.misc.column_flip(data, data['flip'], 'chromosome_1', 'chromosome_2')
    utils.misc.column_flip(data, data['flip'], 'strand_1', 'strand_2')
    utils.misc.column_flip(data, data['flip'], 'position_1', 'position_2')

    # Calculate advancement due to homology based on strand
    data['advance_1'] = np.where(data['strand_1'] == "+", data['homology'], -data['homology'])
    data['advance_2'] = np.where(data['strand_2'] == "+", data['homology'], -data['homology'])

    # Positions have already been flipped, thus we have to
    # add to the second and subtract from the first
    data.loc[data['flip'], 'position_1'] -= data['advance_1']
    data.loc[data['flip'], 'position_2'] += data['advance_2']

    # Track which end of the cluster was the seed end for this read
    data['seed_end'] = np.where(data['flip'], 1-data['read_end'], data['read_end'])

    def revcomp_inserted(row):
        if row['seed_end'] == 0:
            return row['inserted']
        else:
            return utils.misc.reverse_complement(row['inserted'])

    data['inserted'] = data.apply(revcomp_inserted, axis=1)

    data['inslen'] = data['inserted'].apply(len)

    data.set_index(['cluster_id', 'position_1', 'position_2', 'homology', 'inslen'], inplace=True)
    data = data.sort_index()

    agg_f = {'score':sum, 'read_id':len,
             'chromosome_1':max, 'chromosome_2':max,
             'strand_1':max, 'strand_2':max}
    split_data = data.groupby(level=[0, 1, 2, 3, 4])\
                     .agg(agg_f)\
                     .rename(columns={'read_id':'count'})

    selected = split_data['score'].groupby(level=[0]).idxmax().values

    split_data = split_data.loc[selected]

    def calculate_consensus(dat):
        inserted = np.array([np.array(list(a)) for a in dat.values])
        consensus = list()
        for nt_list in inserted.T:
            consensus.append(collections.Counter(nt_list).most_common(1)[0][0])
        consensus = ''.join(consensus)
        return consensus

    split_data['inserted'] = data.loc[selected, 'inserted']\
                                 .groupby(level=[0, 1, 2, 3, 4])\
                                 .apply(calculate_consensus)

    split_data.reset_index(inplace=True)

    return split_data


def predict_breaks(clusters_filename, spanning_filename, split_filename, breakpoints_filename):

    clusters = pd.read_csv(clusters_filename, sep='\t', names=cluster_fields)

    if len(clusters.index) == 0:
        with open(breakpoints_filename, 'w'):
            pass
        return

    split = pd.read_csv(split_filename, sep='\t', names=split_fields, na_values=['.'],
                        converters={'chromosome_1':str, 'chromosome_2':str})
    split['inserted'] = split['inserted'].fillna('')

    spanning = pd.read_csv(spanning_filename, sep='\t', names=spanning_fields,
                           converters={'chromosome':str})

    predictions_0 = predict_breaks_spanning(clusters, spanning)
    predictions_0['prediction_id'] = 0

    predictions_1 = predict_breaks_split(clusters, split)
    predictions_1['prediction_id'] = 1

    predictions = pd.concat([predictions_0, predictions_1], ignore_index=True)
    predictions.sort('cluster_id', inplace=True)

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

    breakpoints = pd.read_csv(breakpoints_filename, sep='\t', names=breakpoint_fields,
                              converters={'chromosome_1':str, 'chromosome_2':str},
                              na_values=['.'])

    breakpoints['inserted'] = breakpoints['inserted'].fillna('')

    breakpoints['inslen'] = breakpoints['inserted'].apply(len)

    data = pd.read_csv(realignments_filename, sep='\t', names=realignment_fields)

    assert data.duplicated(subset=['cluster_id', 'prediction_id', 'library_id', 'read_id', 'read_end']).sum() == 0

    data = data.merge(score_stats, on='aligned_length')

    # Alignment score likelihood and CDF
    data['max_score'] = match_score * data['aligned_length']
    data['score_diff'] = data['max_score'] - data['score']
    data['score_log_likelihood'] = np.log(data['expon_lda']) - \
                                          data['expon_lda'] * data['score_diff']
    data['score_log_cdf'] = -data['expon_lda'] * data['score_diff']

    # Unstack on cluster end
    index_fields = ['cluster_id', 'prediction_id', 'library_id', 'read_id']
    unstack_field = ['cluster_end']
    data_fields = ['read_end', 'aligned_length', 'template_length',
                   'score', 'score_log_likelihood', 'score_log_cdf']
    data = data.set_index(index_fields + unstack_field)[data_fields].unstack()
    data.columns = ['_'.join((a, str(b+1))) for a, b in data.columns.values]
    data.reset_index(inplace=True)

    # Merge insert length from breakpoint predictions
    data = data.merge(breakpoints[['cluster_id', 'prediction_id', 'inslen']],
                      on=['cluster_id', 'prediction_id'])

    data['template_length'] = data['template_length_1'] + data['template_length_2'] + data['inslen']

    # Template length likelihood and CDF
    constant = 1. / ((2 * np.pi)**0.5 * fragment_stddev)
    data['length_z_score'] = (data['template_length'] - fragment_mean) / fragment_stddev
    data['length_log_likelihood'] = -np.log(constant) - np.square(data['length_z_score']) / 2.
    data['length_log_cdf'] = np.log(2. * scipy.stats.norm.sf(data['length_z_score'].abs()))

    data['log_likelihood'] = data['score_log_likelihood_1'] + \
                             data['score_log_likelihood_2'] + \
                             data['length_log_likelihood']
    data['log_cdf'] = data['score_log_cdf_1'] + \
                      data['score_log_cdf_2'] + \
                      data['length_log_cdf']

    data = data.reset_index()

    data = data[likelihoods_fields]

    data.to_csv(likelihoods_filename, sep='\t', index=False, header=False)


def read_merge_write(in_filename, in_names, to_merge, merge_cols, out_filename):
    csv_iter = pd.read_csv(in_filename, sep='\t', names=in_names,
                           iterator=True, chunksize=1000)
    first = True
    for chunk in csv_iter:
        chunk = chunk.merge(to_merge, left_on=merge_cols,
                            right_index=True, how='inner')
        chunk.to_csv(out_filename, sep='\t', mode=('a', 'w')[first], index=False, header=False)
        first = False


def select_clusters(clusters_filename,
                    breakpoints_filename, selected_breakpoints_filename,
                    likelihoods_filename, selected_likelihoods_filename):

    clusters = pd.read_csv(clusters_filename, sep='\t', names=cluster_fields,
                           usecols=['cluster_id', 'library_id', 'read_id'])

    likelihoods_index = ['cluster_id', 'library_id', 'read_id']
    selected = clusters[likelihoods_index].drop_duplicates().set_index(likelihoods_index)

    read_merge_write(likelihoods_filename, likelihoods_fields, selected,
                     ['cluster_id', 'library_id', 'read_id'],
                     selected_likelihoods_filename)

    breakpoints_index = ['cluster_id']
    selected = clusters[breakpoints_index].drop_duplicates().set_index(breakpoints_index)

    read_merge_write(breakpoints_filename, breakpoint_fields, selected,
                     ['cluster_id'],
                     selected_breakpoints_filename)


def select_predictions(breakpoints_filename, selected_breakpoints_filename,
                       likelihoods_filename, selected_likelihoods_filename):

    likelihoods = pd.read_csv(likelihoods_filename, sep='\t', names=likelihoods_fields,
                              usecols=['cluster_id', 'prediction_id', 'log_likelihood'])

    selected = likelihoods.set_index(['cluster_id', 'prediction_id'])\
                          .groupby(level=[0])['log_likelihood']\
                          .idxmax()
    selected = pd.DataFrame(list(selected.values), columns=['cluster_id', 'prediction_id'])

    selected.set_index(['cluster_id', 'prediction_id'], inplace=True)

    read_merge_write(likelihoods_filename, likelihoods_fields, selected,
                     ['cluster_id', 'prediction_id'],
                     selected_likelihoods_filename)

    read_merge_write(breakpoints_filename, breakpoint_fields, selected,
                     ['cluster_id', 'prediction_id'],
                     selected_breakpoints_filename)



