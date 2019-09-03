import sys
import tarfile
import numpy as np
import pandas as pd

import destruct.utils.plots


def load_align_data_hist(filename):

    data = pd.read_table(filename, sep='\t', names=['aligned_length', 'score'])

    align_data_hist = data.groupby(['aligned_length', 'score']).size()
    align_data_hist.name = 'freq'

    return align_data_hist.reset_index()


def create_score_stats(true_scores_filename, match_score, score_stats_filename):
    ''' Infer distribution of null alignment scores and true alignment scores from samples of null and true scores.
    The null samples may include some true scores so run a quick EM mixture model on the null samples to identify
    the actual distribution of null scores.  Output distributions for each alignment length.
    '''
    all_true_scores = load_align_data_hist(true_scores_filename)

    all_true_scores['penalty'] = match_score * all_true_scores['aligned_length'] - all_true_scores['score']

    score_stats = list()

    for aligned_length in all_true_scores['aligned_length'].unique():

        true_scores = all_true_scores.loc[all_true_scores['aligned_length'] == aligned_length]

        expon_lda = 1.0 / np.average(true_scores['penalty'], weights=true_scores['freq'])

        expon_lda = min(expon_lda, 1.0)

        score_stats.append((aligned_length, expon_lda))

    score_stats = pd.DataFrame(score_stats, columns=['aligned_length', 'expon_lda'])

    score_stats.to_csv(score_stats_filename, sep='\t', header=False, index=False)

