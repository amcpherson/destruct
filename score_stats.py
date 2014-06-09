import sys
import tarfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import utils.plots

def load_align_data_hist(filename):
    reader = pd.read_table(filename, sep='\t', names=['aligned_length', 'score'], chunksize=1024*1024)
    align_data_hist = pd.Series()
    for chunk in reader:
        chunk_hist = chunk.groupby(['aligned_length', 'score']).size()
        align_data_hist, chunk_hist = align_data_hist.align(chunk_hist, fill_value=0)
        align_data_hist += chunk_hist
    align_data_hist.name = 'freq'
    return align_data_hist.reset_index()

def create_score_stats(true_scores_filename, match_score, score_stats_filename, plots_tar_filename, library_id):
    ''' Infer distribution of null alignment scores and true alignment scores from samples of null and true scores.
    The null samples may include some true scores so run a quick EM mixture model on the null samples to identify 
    the actual distribution of null scores.  Output distributions for each alignment length.
    '''

    with tarfile.open(plots_tar_filename, 'w') as plots_tar:

        all_true_scores = load_align_data_hist(true_scores_filename)

        all_true_scores['penalty'] = match_score * all_true_scores['aligned_length'] - all_true_scores['score']

        score_stats = list()

        for aligned_length in all_true_scores['aligned_length'].unique():

            true_scores = all_true_scores.loc[all_true_scores['aligned_length'] == aligned_length]

            expon_lda = 1.0 / np.average(true_scores['penalty'], weights=true_scores['freq'])

            score_stats.append((aligned_length, expon_lda))

            fig = plt.figure(figsize=(8,8))

            xmin = true_scores['penalty'].min()
            xmax = true_scores['penalty'].max()
            plot_scores = np.arange(xmin, xmax+1.0, 1.0)

            expon_pdf = lambda x: expon_lda * np.exp(-expon_lda*x)

            utils.plots.filled_density_weighted(plt.gca(), true_scores['penalty'].values, true_scores['freq'].values, 'blue', 0.5, xmin, xmax, 2.0)
            plt.plot(plot_scores, expon_pdf(plot_scores), 'red')
            plt.title('true/null alignment scores for library {0}, aligned length {1}'.format(library_id, aligned_length))
            utils.plots.savefig_tar(plots_tar, fig, 'score_stats_{0}_aligned_length_{1}.pdf'.format(library_id, aligned_length))
            plt.clf()

        score_stats = pd.DataFrame(score_stats, columns=['aligned_length', 'expon_lda'])

        score_stats.to_csv(score_stats_filename, sep='\t', header=False, index=False)

