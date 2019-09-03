
import collections
import pandas as pd
import numpy as np
import scipy
import scipy.stats

import destruct.utils.misc
import destruct.utils.streaming


cluster_fields = ['cluster_id', 'cluster_end',
                  'library_id', 'read_id', 'read_end', 'align_id']


spanning_fields = ['library_id', 'read_id', 'read_end', 'align_id',
                   'chromosome', 'strand', 'position', 'aligned_length',
                   'mate_length', 'mate_score']


split_fields = ['library_id', 'read_id', 'read_end',
                'align_id_1', 'chromosome_1', 'strand_1', 'position_1',
                'align_id_2', 'chromosome_2', 'strand_2', 'position_2',
                'homology', 'inserted', 'score']


breakpoint_fields = ['cluster_id', 'breakpoint_id',
                     'chromosome_1', 'strand_1', 'position_1',
                     'chromosome_2', 'strand_2', 'position_2',
                     'homology', 'count', 'inserted',
                     'mate_score']


realignment_fields = ['cluster_id', 'breakpoint_id', 'cluster_end',
                      'library_id', 'read_id', 'read_end', 'align_id',
                      'aligned_length', 'template_length', 'score']


score_stats_fields = ['aligned_length', 'expon_lda']


likelihoods_fields = ['cluster_id', 'breakpoint_id',
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


def calculate_mate_score(clusters, spanning):

    # Create a table of spanning alignments for these clusters
    data = spanning.merge(clusters, on=['library_id', 'read_id', 'read_end', 'align_id'])

    data = data.groupby(['cluster_id', 'library_id', 'read_id'])[['mate_score']].max()\
               .groupby(level=0).mean()\
               .reset_index()

    return data


def predict_breaks_spanning(clusters, spanning):

    # Create a table of spanning alignments for these clusters
    data = spanning.merge(clusters, on=['library_id', 'read_id', 'read_end', 'align_id'])

    # Add start and end based on position and aligned length
    data['start'] = np.where(data['strand'] == '+', data['position'], data['position'] - data['aligned_length'] + 1)
    data['end'] = np.where(data['strand'] == '+', data['position'] + data['aligned_length'] - 1, data['position'])

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
    data['breakpoint_id'] = 0

    return data


def predict_breaks_split(clusters, split, max_predictions_per_cluster=10):

    # Unstack clusters file by read end, in preparation for
    # merge with split alignments
    index_fields = ['cluster_id', 'library_id', 'read_id']
    data_fields = ['cluster_end', 'align_id']
    paired = pd.merge(clusters.loc[clusters['read_end'] == 0, index_fields + data_fields].drop_duplicates(),
                      clusters.loc[clusters['read_end'] == 1, index_fields + data_fields].drop_duplicates(),
                      on=index_fields,
                      suffixes=('_1', '_2'))

    # Keep track of which reads have ends flipped relative to clusters
    paired['flip'] = paired['cluster_end_1'] != 0
    paired = paired.drop(['cluster_end_1', 'cluster_end_2'], axis=1)
    paired = paired.reset_index()

    # Merge split read information into clusters
    data = split.merge(paired, on=['library_id', 'read_id', 'align_id_1', 'align_id_2'])
    data['position_1'] = data['position_1'].astype(int)
    data['position_2'] = data['position_2'].astype(int)

    # Check for empty split alignments
    if len(data.index) == 0:
      return pd.DataFrame()

    # Flip columns to make them consistent across clusters
    destruct.utils.misc.column_flip(data, data['flip'], 'chromosome_1', 'chromosome_2')
    destruct.utils.misc.column_flip(data, data['flip'], 'strand_1', 'strand_2')
    destruct.utils.misc.column_flip(data, data['flip'], 'position_1', 'position_2')

    # Track which end of the cluster was the seed end for this read
    data['seed_end'] = np.where(data['flip'], 1-data['read_end'], data['read_end'])

    def revcomp_inserted(row):
        if row['seed_end'] == 0:
            return row['inserted']
        else:
            return destruct.utils.misc.reverse_complement(row['inserted'])

    data['inserted'] = data.apply(revcomp_inserted, axis=1)

    data['inslen'] = data['inserted'].apply(len)

    data.set_index(['cluster_id', 'position_1', 'position_2', 'homology', 'inslen'], inplace=True)

    def calculate_consensus(dat):
        inserted = np.array([np.array(list(a)) for a in dat.values])
        consensus = list()
        for nt_list in inserted.T:
            consensus.append(collections.Counter(nt_list).most_common(1)[0][0])
        consensus = ''.join(consensus)
        return consensus

    agg_f = {'score':sum, 'read_id':len,
             'chromosome_1':max, 'chromosome_2':max,
             'strand_1':max, 'strand_2':max,
             'inserted':calculate_consensus}
    split_data = data.groupby(level=[0, 1, 2, 3, 4])\
                     .agg(agg_f)\
                     .rename(columns={'read_id':'count'})\
                     .reset_index()

    split_data = split_data.sort_values(['cluster_id', 'score'])

    split_data['breakpoint_id'] = split_data.groupby('cluster_id').cumcount(ascending=False)

    split_data = split_data[split_data['breakpoint_id'] < max_predictions_per_cluster]

    split_data['breakpoint_id'] += 1

    return split_data


def predict_breaks(clusters_filename, spanning_filename, split_filename, breakpoints_filename):

    clusters = pd.read_csv(clusters_filename, sep='\t', names=cluster_fields)

    if len(clusters.index) == 0:
        with open(breakpoints_filename, 'w'):
            pass
        return

    merge_columns = ['library_id', 'read_id', 'read_end', 'align_id']

    clusters_alignments = clusters[merge_columns].drop_duplicates()

    split_iter = pd.read_csv(split_filename, sep='\t', names=split_fields,
                             converters={'chromosome_1':str, 'chromosome_2':str, 'inserted':str},
                             iterator=True, chunksize=1000000)

    split_merge_columns = {'1':['library_id', 'read_id', 'read_end', 'align_id_1'],
                           '2':['library_id', 'read_id', 'read_end', 'align_id_2']}

    def filter_split(df):
        filtered = list()
        for side, left_merge_columns in split_merge_columns.items():
            df_2 = pd.merge(df, clusters_alignments,
                            left_on=left_merge_columns,
                            right_on=merge_columns,
                            how='inner')
            df_2 = df_2.drop(['align_id'], axis=1)
            filtered.append(df_2)
        return pd.concat(filtered).drop_duplicates()

    split = pd.concat([filter_split(chunk) for chunk in split_iter])

    split.loc[split['inserted'] == '.', 'inserted'] = ''

    spanning_iter = pd.read_csv(spanning_filename, sep='\t', names=spanning_fields,
                                converters={'chromosome':str},
                                iterator=True, chunksize=1000000)

    def filter_spanning(df):
        return pd.merge(df, clusters_alignments, how='inner')

    spanning = pd.concat([filter_spanning(chunk) for chunk in spanning_iter])

    predictions_0 = predict_breaks_spanning(clusters, spanning)

    predictions_1 = predict_breaks_split(clusters, split)

    predictions = pd.concat([predictions_0, predictions_1], ignore_index=True)
    predictions.sort_values('cluster_id', inplace=True)

    mate_scores = calculate_mate_score(clusters, spanning)
    predictions = predictions.merge(mate_scores, on='cluster_id')

    predictions.loc[predictions['inserted'] == '', 'inserted'] = '.'

    predictions = predictions[breakpoint_fields]
    predictions.to_csv(breakpoints_filename, sep='\t', index=False, header=False)


def calculate_cluster_weights(breakpoints_filename, weights_filename):

    epsilon = 0.0001
    itx_distance = 1000000000

    breakpoints = pd.read_csv(breakpoints_filename, sep='\t', names=breakpoint_fields,
                              converters={'chromosome_1':str, 'chromosome_2':str, 'inserted':str})

    breakpoints = breakpoints[breakpoints['breakpoint_id'] == 0]

    breakpoints['distance'] = np.absolute(breakpoints['position_1'] - breakpoints['position_2']) + 1.0
    breakpoints.loc[breakpoints['chromosome_1'] != breakpoints['chromosome_2'], 'distance'] = itx_distance
    breakpoints['weight'] = 1.0 + epsilon * np.log(breakpoints['distance'])

    breakpoints = breakpoints.sort_values('cluster_id')
    breakpoints[['cluster_id', 'weight']].to_csv(weights_filename, sep='\t', index=False, header=False)


def calculate_realignment_likelihoods(breakpoints_filename, realignments_filename, score_stats_filename,
                                      likelihoods_filename, match_score, fragment_mean, fragment_stddev):

    match_score = float(match_score)
    fragment_mean = float(fragment_mean)
    fragment_stddev = float(fragment_stddev)

    score_stats = pd.read_csv(score_stats_filename, sep='\t', names=score_stats_fields)

    breakpoints = pd.read_csv(breakpoints_filename, sep='\t', names=breakpoint_fields,
                              converters={'chromosome_1':str, 'chromosome_2':str, 'inserted':str})

    breakpoints.loc[breakpoints['inserted'] == '.', 'inserted'] = ''

    breakpoints['inslen'] = breakpoints['inserted'].apply(len)

    data = pd.read_csv(realignments_filename, sep='\t', names=realignment_fields)

    if len(data.index) == 0:
        with open(likelihoods_filename, 'w'):
            pass
        return

    data = data.merge(score_stats, on='aligned_length')

    # Alignment score likelihood and CDF
    data['max_score'] = match_score * data['aligned_length']
    data['score_diff'] = data['max_score'] - data['score']
    data['score_log_likelihood'] = np.log(data['expon_lda']) - \
                                          data['expon_lda'] * data['score_diff']
    data['score_log_cdf'] = -data['expon_lda'] * data['score_diff']

    # Unstack on cluster end
    index_fields = ['cluster_id', 'breakpoint_id', 'library_id', 'read_id']
    data_fields = ['read_end', 'aligned_length', 'template_length',
                   'score', 'score_log_likelihood', 'score_log_cdf']
    data = pd.merge(data.loc[data['cluster_end'] == 0, index_fields + data_fields].drop_duplicates(),
                    data.loc[data['cluster_end'] == 1, index_fields + data_fields].drop_duplicates(),
                    on=index_fields,
                    suffixes=('_1', '_2'))

    # For foldbacks and similar, a read could conceivable align to both sides of the breakpoint
    # remove situations in which one read end is aligned to both cluster ends
    data = data[data['read_end_1'] != data['read_end_2']]

    # Merge insert length from breakpoint predictions
    data = data.merge(breakpoints[['cluster_id', 'breakpoint_id', 'inslen']],
                      on=['cluster_id', 'breakpoint_id'])

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

    index_fields = ['cluster_id', 'breakpoint_id', 'library_id', 'read_id']
    data = data.sort_values(index_fields + ['log_likelihood'])\
               .groupby(index_fields)\
               .last()\
               .reset_index()

    data = data[likelihoods_fields]

    data.to_csv(likelihoods_filename, sep='\t', index=False, header=False)


def select_clusters(clusters_filename,
                    breakpoints_filename, selected_breakpoints_filename,
                    likelihoods_filename, selected_likelihoods_filename):

    clusters = pd.read_csv(clusters_filename, sep='\t', names=cluster_fields,
                           usecols=['cluster_id', 'library_id', 'read_id'])
    clusters = clusters.drop_duplicates()

    breakpoints_iter = pd.read_csv(breakpoints_filename, sep='\t', names=breakpoint_fields,
                                   converters={'chromosome_1':str, 'chromosome_2':str, 'inserted':str},
                                   iterator=True, chunksize=1000000)

    cluster_ids = clusters[['cluster_id']].drop_duplicates()

    destruct.utils.streaming.read_select_write(breakpoints_iter, cluster_ids, selected_breakpoints_filename)

    likelihoods_iter = pd.read_csv(likelihoods_filename, sep='\t', names=likelihoods_fields,
                                   iterator=True, chunksize=1000000)

    destruct.utils.streaming.read_select_write(likelihoods_iter, clusters, selected_likelihoods_filename)


def select_breakpoint_prediction(likelihoods, template_length_min_threshold):
    """ Select and filter breakpoint predictions.

    Args:
        likelihoods(pandas.DataFrame): likelihoods table
        template_length_min_threshold (int): min template length filter

    Select the maximum likelihoods breakpoint prediction, preferring solutions
    with split reads.  Filter predictions based on minimum of template lengths on
    either side of the breakpoint.

    """

    # Calculate total breakpoint likelihood, max template lengths
    agg_f = {
        'log_likelihood':sum,
        'template_length_1':max,
        'template_length_2':max,
    }
    data = likelihoods.groupby(['cluster_id', 'breakpoint_id']).agg(agg_f).reset_index()

    # Select highest likelihood breakpoint predictions
    # Prefer a higher breakpoint id, thus preferring solutions with split reads
    data.sort_values(['cluster_id', 'log_likelihood', 'breakpoint_id'],
        ascending=[True, False, False], inplace=True)
    selected = data.groupby('cluster_id', sort=False).first().reset_index()

    # Filter by template length
    selected['template_length_min'] = selected[['template_length_1', 'template_length_2']].min(axis=1)
    selected.drop(['template_length_1', 'template_length_2'], axis=1, inplace=True)
    selected = selected[selected['template_length_min'] >= template_length_min_threshold]

    return selected


def select_predictions(breakpoints_filename, selected_breakpoints_filename,
                       likelihoods_filename, selected_likelihoods_filename,
                       mate_score_threshold, template_length_min_threshold,
                       min_alignment_log_likelihood):

    read_likelihoods_iter = pd.read_csv(likelihoods_filename, sep='\t', names=likelihoods_fields,
        usecols=['cluster_id', 'breakpoint_id', 'log_likelihood', 'template_length_1', 'template_length_2'],
        chunksize=int(1e7))

    read_likelihoods_iter = destruct.utils.streaming.group_aware_iter(read_likelihoods_iter, ['cluster_id'])

    selected = pd.concat([
        select_breakpoint_prediction(df, template_length_min_threshold)
        for df in read_likelihoods_iter], ignore_index=True)

    mate_score = pd.read_csv(breakpoints_filename, sep='\t', names=breakpoint_fields,
        converters={'chromosome_1':str, 'chromosome_2':str, 'inserted':str},
        usecols=['cluster_id', 'breakpoint_id', 'mate_score'])

    mate_score = mate_score.loc[mate_score['mate_score'] <= mate_score_threshold, ['cluster_id', 'breakpoint_id']]
    selected = selected.merge(mate_score)

    selected = selected[['cluster_id', 'breakpoint_id']].drop_duplicates()

    breakpoints_iter = pd.read_csv(breakpoints_filename, sep='\t', names=breakpoint_fields,
                                   converters={'chromosome_1':str, 'chromosome_2':str, 'inserted':str},
                                   iterator=True, chunksize=1000000)

    destruct.utils.streaming.read_select_write(breakpoints_iter, selected, selected_breakpoints_filename)

    likelihoods_iter = pd.read_csv(likelihoods_filename, sep='\t', names=likelihoods_fields,
                                   iterator=True, chunksize=1000000)

    def likelihoods_filter(chunk):
        chunk = chunk.merge(selected, how='inner')
        chunk = chunk[chunk['log_likelihood'] >= min_alignment_log_likelihood]
        return chunk

    destruct.utils.streaming.read_filter_write(likelihoods_iter, likelihoods_filter, selected_likelihoods_filename)


