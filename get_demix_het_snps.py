import os
import sys
import re
import argparse
import pandas as pd


argparser = argparse.ArgumentParser()
argparser.add_argument('demix_tmp', help='deMix temporary directory')
argparser.add_argument('segment', help='Chromosomal segment as chr:start-end')
argparser.add_argument('results', help='Table of haplotyped het snps')
args = argparser.parse_args()


try:
    chromosome, start, end = re.match('(.*):(.*)-(.*)', args.segment).groups()
    start = int(start)
    end = int(end)
except:
    raise Exception('unable to parse segment ' + args.segment)

haps_filename = os.path.join(args.demix_tmp, 'tmp', 'haps.'+chromosome)

try:
    haps = pd.read_csv(haps_filename, sep='\t')
except IOError:
    raise Exception('unable to open ' + haps_filename + ' you may have the wrong temps directory, or temps were cleaned up')

haps = haps[(haps['pos'] >= start) & (haps['pos'] <= end)]

haps.set_index('allele_label', inplace=True)
haps['hap_size'] = haps.groupby(level=0).size()
haps.reset_index(inplace=True)

haps.to_csv(args.results, sep='\t', index=False)

