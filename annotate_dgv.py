import pygenes
from collections import *
import csv
import sys

dgv_filename = sys.argv[1]

class DGVDatabase(object):
    def __init__(self, dgv_filename):
        self.variations = dict()
        chrvars = defaultdict(list)
        with open(dgv_filename, 'r') as dgv_file:
            dgv_reader = csv.reader(dgv_file, delimiter='\t')
            dgv_header = next(dgv_reader)
            for row in dgv_reader:
                assert row[0].startswith('Variation_')
                id = int(row[0][len('Variation_'):])
                chr = row[2][3:]
                start = int(row[3])
                end = int(row[4])
                chrvars[chr].append((id, start, end))
                self.variations[id] = (start, end)
        self.intervals = dict()
        for chr, vars in chrvars.iteritems():
            self.intervals[chr] = pygenes.IntervalTree(vars)
    def query(self, region):
        if region[0][0] != region[1][0]:
            return
        if region[0][0] not in self.intervals:
            return
        break1 = (region[0][2], region[0][3])[region[0][1] == '+']
        break2 = (region[1][2], region[1][3])[region[1][1] == '+']
        start, end = sorted([break1, break2])
        ids = self.intervals[region[0][0]].find_overlapping(start, end)
        for idx in range(0, len(ids)):
            startdiff = abs(start - self.variations[ids[idx]][0])
            enddiff = abs(end - self.variations[ids[idx]][1])
            if startdiff < 500 and enddiff < 500:
                yield 'Variation_' + str(ids[idx])

dgv = DGVDatabase(dgv_filename)

reader = csv.reader(sys.stdin, delimiter='\t')
header = next(reader)
print '\t'.join(header + ['dgv'])
for row in reader:
    fields = dict(zip(header, row))
    region = [[fields['chromosome_1'], fields['strand_1'], int(fields['start_1']), int(fields['end_1'])],
              [fields['chromosome_2'], fields['strand_2'], int(fields['start_2']), int(fields['end_2'])]]
    print '\t'.join(row + [','.join(dgv.query(region))])
