
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import random
import itertools
import seaborn


chromosomes = [str(a) for a in range(1, 23)] + ['X']
chromosome_indices = dict([(chromosome, idx) for idx, chromosome in enumerate(chromosomes)])

cnv = pd.read_csv('/Users/amcphers/Analysis/demix_interactive/patient_3.preds.tsv', sep='\t', converters={'chr':str})

cnv = cnv.loc[(cnv['library_id'] == 'adnexa_site_1')]
cnv = cnv.loc[(cnv['chr'].isin(chromosomes))]

cnv['chr_index'] = cnv['chr'].apply(lambda a: chromosome_indices[a])

cnv = cnv.sort(['chr_index', 'start'])

chromosome_maximum = cnv.groupby('chr_index')['end'].max()

chromosome_end = np.cumsum(chromosome_maximum)
chromosome_start = chromosome_end.shift(1)
chromosome_start[0] = 0

cnv.set_index('chr_index', inplace=True)
cnv['chromosome_start'] = chromosome_start
cnv['chromosome_end'] = chromosome_end
cnv.reset_index(inplace=True)

cnv['chromosome_mid'] = 0.5 * (cnv['chromosome_start'] + cnv['chromosome_end'])

cnv['start'] += cnv['chromosome_start']
cnv['end'] += cnv['chromosome_start']

mingap = 1000

copies_max = 4.0





fig = plt.figure(figsize=(16,16))

gs = matplotlib.gridspec.GridSpec(2, 1, height_ratios=(4, 1))

ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

color_set = plt.get_cmap('Set1')
color_set = [color_set(float(i)/len(chromosomes)) for i in range(len(chromosomes))]
chromosome_color = lambda c: color_set[chromosomes.index(c)]
cs = [chromosome_color(c) for c in cnv['chr'].values]
            
major_minor_scatter = ax1.scatter(cnv['major_raw'], cnv['minor_raw'],
                                  s=cnv['length']/20000.0, 
                                  facecolor=cs, edgecolor=cs, linewidth=0.0,
                                  picker=True)

ax1.set_xlim((-0.5, copies_max))
ax1.set_ylim((-0.5, 0.8*copies_max))

lgnd = ax1.legend([plt.Circle((0, 0), radius=1, color=chromosome_color(c), picker=True) for c in chromosomes], chromosomes, loc=2)
lgnd_patches = list(lgnd.get_patches())

for patch in lgnd_patches:
    patch.set_picker(True)

major_segments = list()
minor_segments = list()
major_connectors = list()
minor_connectors = list()

for (idx, row), (next_idx, next_row) in itertools.izip_longest(cnv.iterrows(), cnv.iloc[1:].iterrows(), fillvalue=(None, None)):
    major_segments.append([(row['start'], row['major_raw']), (row['end'], row['major_raw'])])
    minor_segments.append([(row['start'], row['minor_raw']), (row['end'], row['minor_raw'])])
    if next_row is not None and next_row['start'] - row['end'] < mingap:
        major_connectors.append([(row['end'], row['major_raw']), (next_row['start'], next_row['major_raw'])])
        minor_connectors.append([(row['end'], row['minor_raw']), (next_row['start'], next_row['minor_raw'])])

major_segments = matplotlib.collections.LineCollection(major_segments, colors='r')
minor_segments = matplotlib.collections.LineCollection(minor_segments, colors='b')
major_connectors = matplotlib.collections.LineCollection(major_connectors, colors='r')
minor_connectors = matplotlib.collections.LineCollection(minor_connectors, colors='b')

major_segments.set_picker(True)
minor_segments.set_picker(True)

linewidths = np.array([1] * len(cnv.index))
major_minor_scatter.set_linewidths(linewidths)
major_segments.set_linewidths(linewidths)
minor_segments.set_linewidths(linewidths)

scatter_edgecolors = np.array(['b'] * len(cnv.index))
major_minor_scatter.set_edgecolors(scatter_edgecolors)

ax2.add_collection(major_segments)
ax2.add_collection(minor_segments)
ax2.add_collection(major_connectors)
ax2.add_collection(minor_connectors)
ax2.set_xlim((cnv['start'].min(), cnv['end'].max()))
ax2.set_ylim((-0.2, copies_max + 0.2))

ax2.set_xticks([0] + sorted(cnv['chromosome_end'].unique()))
ax2.set_xticklabels([])

ax2.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(sorted(cnv['chromosome_mid'].unique())))
ax2.xaxis.set_minor_formatter(matplotlib.ticker.FixedFormatter(chromosomes))

ax2.grid(False, which="minor")

def onpick1(event):
    print event.artist
    if isinstance(event.artist, matplotlib.patches.Rectangle):
        for ind, patch in enumerate(lgnd_patches):
            if patch == event.artist:
                print chromosomes[ind]
    elif isinstance(event.artist, matplotlib.collections.PathCollection):
        pass
    elif isinstance(event.artist, matplotlib.collections.LineCollection):
        linewidths = np.array([1] * len(cnv.index))
        linewidths[event.ind] = 4
        scatter_edgecolors = np.array(['b'] * len(cnv.index))
        scatter_edgecolors[event.ind] = 'yellow'
        major_minor_scatter.set_edgecolors(scatter_edgecolors)
        major_minor_scatter.set_linewidths(linewidths)
        major_segments.set_linewidths(linewidths)
        minor_segments.set_linewidths(linewidths)
    event.canvas.draw()
        
fig.canvas.mpl_connect('pick_event', onpick1)


plt.show()

