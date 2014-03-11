import sys
import numpy as np
import pandas as pd
from scipy.stats import nbinom

args = sys.argv[1:]

def load_align_data_hist(filename):
    reader = pd.read_table(filename, sep='\t', names=['aligned_length', 'score'], chunksize=1024*1024)
    align_data_hist = pd.Series()
    for chunk in reader:
        chunk_hist = chunk.groupby(['aligned_length', 'score']).size()
        align_data_hist, chunk_hist = align_data_hist.align(chunk_hist, fill_value=0)
        align_data_hist += chunk_hist
    align_data_hist.name = 'freq'
    return align_data_hist.reset_index()

all_true_scores = load_align_data_hist(args[0])
all_null_scores = load_align_data_hist(args[1])
match_score = int(args[2])

all_true_scores['penalty'] = match_score * all_true_scores['aligned_length'] - all_true_scores['score']
all_null_scores['penalty'] = match_score * all_null_scores['aligned_length'] - all_null_scores['score']

score_stats = list()

for aligned_length in all_true_scores['aligned_length'].unique():

    true_scores = all_true_scores.loc[all_true_scores['aligned_length'] == aligned_length]
    null_scores = all_null_scores.loc[all_null_scores['aligned_length'] == aligned_length]

    mean_zd1 = np.average(true_scores['penalty'], weights=true_scores['freq'])
    var_zd1 = np.average((true_scores['penalty']-mean_zd1)**2, weights=true_scores['freq'])
    size_zd1 = mean_zd1 * mean_zd1 / (var_zd1 - mean_zd1)
    prob_zd1 = size_zd1 / (size_zd1 + mean_zd1)

    # Initialize membership
    pdi1 = np.array([0.5]*len(null_scores.index))

    # Initialize pi
    pi = sum(pdi1 * null_scores['freq']) / sum(null_scores['freq'])

    for i in xrange(100):

        mean_zd0 = sum(null_scores['penalty'] * null_scores['freq'] * (1 - pdi1)) / sum(null_scores['freq'] * (1 - pdi1))
        var_zd0 = sum(null_scores['freq'] * (1 - pdi1) * (null_scores['penalty'] - mean_zd0)**2) / sum(null_scores['freq'] * (1 - pdi1))
        size_zd0 = mean_zd0 * mean_zd0 / (var_zd0 - mean_zd0)
        prob_zd0 = size_zd0 / (size_zd0 + mean_zd0)
    
        pdi1 = nbinom.pmf(null_scores['penalty'], size_zd1, prob_zd1) * pi / (nbinom.pmf(null_scores['penalty'], size_zd1, prob_zd1) * pi + nbinom.pmf(null_scores['penalty'], size_zd0, prob_zd0) * (1 - pi))
    
        pi = sum(pdi1 * null_scores['freq']) / sum(null_scores['freq'])
    
    if np.isnan(pi):
        pi = 0.5
        size_zd0 = size_zd1
        prob_zd0 = prob_zd1

    score_stats.append((aligned_length, pi, size_zd1, prob_zd1, size_zd0, prob_zd0))

score_stats = pd.DataFrame(score_stats, columns=['aligned_length', 'prob_true', 'score_true_size', 'score_true_prob', 'score_null_size', 'score_null_prob'])

score_stats.to_csv(args[3], sep='\t', header=False, index=False)

