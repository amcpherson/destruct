import collections
import csv
import errno
import itertools
import os
import tarfile
import gzip
import numpy as np
import pandas as pd
import pypeliner
import pygenes

import destruct.utils.plots
import destruct.utils.seq
import destruct.predict_breaks


def prepare_seed_fastq(reads_1_fastq, reads_2_fastq, seed_length, seed_fastq):
    opener_1 = (open, gzip.open)[reads_1_fastq.endswith('.gz')]
    opener_2 = (open, gzip.open)[reads_2_fastq.endswith('.gz')]
    with opener_1(reads_1_fastq, 'rt') as reads_1, opener_2(reads_2_fastq, 'rt') as reads_2, open(seed_fastq, 'wt') as seed:
        fastq_lines = [[], []]
        for fastq_1_line, fastq_2_line in zip(reads_1, reads_2):
            fastq_lines[0].append(fastq_1_line.rstrip())
            fastq_lines[1].append(fastq_2_line.rstrip())
            if len(fastq_lines[0]) == 4:
                for read_end in (0, 1):
                    if len(fastq_lines[read_end][1]) > seed_length:
                        fastq_lines[read_end][1] = fastq_lines[read_end][1][0:seed_length]
                        fastq_lines[read_end][3] = fastq_lines[read_end][3][0:seed_length]
                    for line in fastq_lines[read_end]:
                        seed.write(line + '\n')
                fastq_lines = [[], []]


class ConcordantReadStats(object):
    def __init__(self, stats, fragment_length_num_stddevs):
        self.stats = stats
        self.fragment_length_num_stddevs = fragment_length_num_stddevs

    @property
    def fragment_length_mean(self):
        return float(self.stats['fragment_mean'])

    @property
    def fragment_length_stddev(self):
        return float(self.stats['fragment_stddev'])

    @property
    def fragment_length_min(self):
        return int(self.fragment_length_mean - self.fragment_length_num_stddevs * self.fragment_length_stddev)

    @property
    def fragment_length_max(self):
        return int(self.fragment_length_mean + self.fragment_length_num_stddevs * self.fragment_length_stddev)


def read_stats(stats_filename, fragment_length_num_stddevs):
    stats = pd.read_csv(stats_filename, sep='\t')
    flen_stats = stats.loc[stats['type'] == 'fragment_length'].drop('type', axis=1)
    flen_stats = flen_stats.astype(float)
    fragment_count = flen_stats['value'].sum()
    fragment_mean = (flen_stats['key'] * flen_stats['value']).sum() / fragment_count
    fragment_variance = ((flen_stats['key'] - fragment_mean) * (flen_stats['key'] - fragment_mean) * flen_stats['value']).sum() / (fragment_count - 1)
    fragment_stddev = fragment_variance**0.5
    return ConcordantReadStats({'fragment_mean': fragment_mean, 'fragment_stddev': fragment_stddev}, fragment_length_num_stddevs)


def write_stats_table(library_ids, lib_stats, stats_table_filename):
    with open(stats_table_filename, 'wt') as stats_table_file:
        for lib_name, library_id in library_ids.items():
            stats_table_file.write(str(library_id) + '\t')
            stats_table_file.write(str(lib_stats[lib_name].fragment_length_mean) + '\t')
            stats_table_file.write(str(lib_stats[lib_name].fragment_length_stddev) + '\n')


def split_file_byline(in_filename, lines_per_file, out_filename_callback):
    with open(in_filename, 'rt') as in_file:
        file_number = 0
        out_file = None
        out_file_lines = None
        try:
            for line in in_file:
                if out_file is None or out_file_lines == lines_per_file:
                    if out_file is not None:
                        out_file.close()
                    out_file = open(out_filename_callback(file_number), 'wt')
                    out_file_lines = 0
                    file_number += 1
                out_file.write(line)
                out_file_lines += 1
        finally:
            if out_file is not None:
                out_file.close()


def split_fastq(in_filename, num_reads_per_file, out_filename_callback):
    with gzip.open(in_filename, 'rt') as in_file:
        file_number = 0
        out_file = None
        out_file_read_count = None
        try:
            for name, seq, comment, qual in itertools.zip_longest(*[in_file]*4):
                if out_file is None or out_file_read_count == num_reads_per_file:
                    if out_file is not None:
                        out_file.close()
                    out_file = open(out_filename_callback(file_number), 'wt')
                    out_file_read_count = 0
                    file_number += 1
                out_file.write(name)
                out_file.write(seq)
                out_file.write(comment)
                out_file.write(qual)
                out_file_read_count += 1
        finally:
            if out_file is not None:
                out_file.close()


def merge_files_by_line(in_filenames, out_filename):
    with open(out_filename, 'wt') as out_file:
        for id, in_filename in sorted(in_filenames.items()):
            with open(in_filename, 'rt') as in_file:
                for line in in_file:
                    out_file.write(line)


def create_library_ids(library_names):
    return dict([(library_name, library_id) for library_id, library_name in enumerate(library_names)])


def merge_alignment_files(in_filenames, out_filename, library_idxs):
    with open(out_filename, 'wt') as out_file:
        for lib_id, in_filename in in_filenames.items():
            idx = library_idxs[lib_id]
            with open(in_filename, 'rt') as in_file:
                for line in in_file:
                    line = str(idx) + line[line.index('\t'):]
                    out_file.write(line)


def merge_sorted_files_by_line(in_filenames, out_filename, temp_space, sort_fields):
    try:
        os.makedirs(temp_space)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    pypeliner.commandline.execute(*['sort', '-T', temp_space, '-m', '-n', '-k', sort_fields] + list(in_filenames.values()) + ['>', out_filename])


def generate_chromosome_args(chromosomes):
    args = list()
    for chromosome_pair in itertools.combinations_with_replacement(chromosomes, 2):
        args.append('--inclchrompair ' + ','.join(chromosome_pair))
    args.append('--exclchrompairs ' + ','.join(chromosomes))
    return dict(enumerate(args))


def read_clusters_breakpoints(clusters_filename, breakpoints_filename):
    with open(clusters_filename, 'rt') as clusters_file, open(breakpoints_filename, 'rt') as breakpoints_file:
        clusters_reader = csv.reader(clusters_file, delimiter='\t')
        breakpoints_reader = csv.reader(breakpoints_file, delimiter='\t')
        cluster_iter = itertools.groupby(clusters_reader, lambda row: row[0])
        breakend_iter = itertools.groupby(breakpoints_reader, lambda row: row[0])
        for (cluster_id_1, cluster_rows), (cluster_id_2, breakend_rows) in zip(cluster_iter, breakend_iter):
            if cluster_id_1 != cluster_id_2:
                raise ValueError('Consistency issue between clusters and breakpoints for ' + clusters_filename + ' and ' + breakpoints_filename)
            yield cluster_id_1, cluster_rows, breakend_rows


def merge_clusters(in_clusters_filenames, in_breakpoints_filenames,
                   out_clusters_filename, out_breakpoints_filename, debug_filename):
    new_cluster_id = 0
    with open(out_clusters_filename, 'wt') as out_clusters_file, \
         open(out_breakpoints_filename, 'wt') as out_breakpoints_file, \
         open(debug_filename, 'wt') as debug_file:
        for idx, in_clusters_filename in in_clusters_filenames.items():
            in_breakpoints_filename = in_breakpoints_filenames[idx]
            for cluster_id, cluster_rows, breakend_rows in read_clusters_breakpoints(in_clusters_filename, in_breakpoints_filename):
                for row in cluster_rows:
                    row[0] = str(new_cluster_id)
                    out_clusters_file.write('\t'.join(row) + '\n')
                for row in breakend_rows:
                    row[0] = str(new_cluster_id)
                    out_breakpoints_file.write('\t'.join(row) + '\n')
                debug_file.write('{0}\t{1}\t{2}\n'.format(new_cluster_id, idx, cluster_id))
                new_cluster_id += 1


def tabulate_reads(clusters_filename, likelihoods_filename, library_ids, reads1_filenames, reads2_filenames, reads_table_filename):
    fields = ['cluster_id', 'cluster_end', 'lib_id', 'read_id', 'read_end', 'align_id']
    clusters = pd.read_csv(clusters_filename, sep='\t', names=fields, usecols=['cluster_id', 'lib_id', 'read_id'])
    clusters = clusters.drop_duplicates().set_index(['lib_id', 'read_id']).sort_index()['cluster_id']

    # The likelihoods file contains a list of read alignments that have passed
    # filtering, use this to generate a list of library / read id pairs that
    # can be used to annotate each read cluster assignment as filtered or not
    likelihoods = pd.read_csv(likelihoods_filename, sep='\t',
                              names=destruct.predict_breaks.likelihoods_fields,
                              converters=converters)
    passed_reads = set(zip(likelihoods.library_id, likelihoods.read_id))

    with open(reads_table_filename, 'wt') as reads_table_file:
        for lib_name in set(reads1_filenames.keys()).union(set(reads2_filenames.keys())):
            lib_id = library_ids[lib_name]
            for reads_filename in [reads1_filenames[lib_name], reads2_filenames[lib_name]]:
                with gzip.open(reads_filename, 'rt') as reads_file:
                    for name, seq, comment, qual in itertools.zip_longest(*[(a.rstrip() for a in reads_file)]*4):
                        assert name[0] == '@'
                        assert name[-1] == '1' or name[-1] == '2'
                        assert name[-2] == '/'
                        fragment_id = int(name[1:-2])
                        read_end = name[-1]
                        if (lib_id, fragment_id) in clusters.index:
                            cluster_id = clusters.loc[(lib_id, fragment_id)]
                            assert np.issubdtype(type(cluster_id), np.integer)
                            filtered = 'False' if (lib_id, fragment_id) in passed_reads else 'True'
                            reads_table_file.write('\t'.join([str(cluster_id), str(lib_id), str(fragment_id), read_end, seq, qual, comment, filtered]) + '\n')


class DGVDatabase(object):
    def __init__(self, dgv_filename):
        self.variations = list()
        chrvars = collections.defaultdict(list)
        with open(dgv_filename, 'rt') as dgv_file:
            dgv_reader = csv.reader(dgv_file, delimiter='\t')
            dgv_header = next(dgv_reader)
            for row in dgv_reader:
                id = row[0]
                chr = row[1]
                start = int(row[2])
                end = int(row[3])
                idx = len(self.variations)
                chrvars[chr].append((idx, start, end))
                self.variations.append((id, start, end))
        self.intervals = dict()
        for chr, vars in chrvars.items():
            self.intervals[chr] = pygenes.IntervalTree(vars)

    def query(self, chromosome, start, end):
        if chromosome not in self.intervals:
            return
        idxs = self.intervals[chromosome].find_overlapping(start, end)
        for idx in [idxs[a] for a in range(0, len(idxs))]:
            startdiff = abs(start - self.variations[idx][1])
            enddiff = abs(end - self.variations[idx][2])
            if startdiff < 500 and enddiff < 500:
                yield self.variations[idx][0]


def merge_tars(output_filename, *input_filename_sets):
    with tarfile.open(output_filename, 'wt') as output_tar:
        for input_filenames in input_filename_sets:
            for input_filename in input_filenames.values():
                with tarfile.open(input_filename, 'rt') as in_tar:
                    for tarinfo in in_tar:
                        output_tar.addfile(tarinfo, in_tar.extractfile(tarinfo))


converters = {'chromosome': str,
              'chromosome_1': str,
              'chromosome_2': str,
              'inserted': str}


def create_sequence(row, reference_sequences):
    breakend_sequences = ['', '']
    expected_strands = ('+', '-')
    inserted = ''
    if inserted != '.':
        inserted = row['inserted']
    for side in (0, 1):
        chromosome = row['chromosome_{0}'.format(side+1)]
        strand = row['strand_{0}'.format(side+1)]
        position = row['position_{0}'.format(side+1)]
        length = row['template_length_{0}'.format(side+1)]
        if strand == '+':
            start = position - length + 1
            end = position
        else:
            start = position
            end = position + length - 1
        breakend_sequences[side] = reference_sequences[chromosome][start-1:end]
        if strand != expected_strands[side]:
            breakend_sequences[side] = destruct.utils.misc.reverse_complement(breakend_sequences[side])
    return breakend_sequences[0] + '[' + inserted + ']' + breakend_sequences[1]


def annotate_genes(row, gene_models):

    for side in (0, 1):

        chromosome = row['chromosome_{0}'.format(side+1)]
        position = row['position_{0}'.format(side+1)]

        nearest_gene_ids = gene_models.find_nearest_genes(chromosome, position)
        gene_id = 'NA'
        gene_name = 'NA'
        gene_location = 'NA'
        if len(nearest_gene_ids) > 0:
            gene_id = nearest_gene_ids[0]
            gene_name = gene_models.get_gene(gene_id).name
            gene_location = gene_models.calculate_gene_location(gene_id, position)

        row['gene_id_{0}'.format(side+1)] = gene_id
        row['gene_name_{0}'.format(side+1)] = gene_name
        row['gene_location_{0}'.format(side+1)] = gene_location

    return row


def query_dgv(row, dgv):

    if row['chromosome_1'] != row['chromosome_2']:
        return 'NA'

    chromosome = row['chromosome_1']
    start, end = sorted((row['position_1'], row['position_2']))

    variants = list(dgv.query(chromosome, start, end))

    if len(variants) == 0:
        return 'NA'

    return ', '.join(variants)


def tabulate_results(breakpoints_filename, likelihoods_filename, library_ids,
                     genome_fasta, gtf_filename, dgv_filename,
                     breakpoint_table, breakpoint_library_table):

    lib_names = pd.DataFrame(library_ids.items(), columns=['library', 'library_id'])

    breakpoints = pd.read_csv(breakpoints_filename, sep='\t',
                              names=destruct.predict_breaks.breakpoint_fields,
                              converters=converters)
    breakpoints = breakpoints.drop(['breakpoint_id'], axis=1)
    breakpoints = breakpoints.rename(columns={'count':'num_split'})
    breakpoints.loc[breakpoints['inserted'] == '.', 'inserted'] = ''

    likelihoods = pd.read_csv(likelihoods_filename, sep='\t',
                              names=destruct.predict_breaks.likelihoods_fields,
                              converters=converters)
    likelihoods = likelihoods.drop(['breakpoint_id'], axis=1)

    breakpoint_reads = (
        likelihoods.groupby(['cluster_id', 'library_id'])
        .size()
        .reset_index()
    )
    breakpoint_reads.columns = ['cluster_id', 'library_id', 'num_reads']

    breakpoint_unique_reads = (
        likelihoods.drop_duplicates(['cluster_id', 'library_id', 'template_length_1', 'template_length_2'])
        .groupby(['cluster_id', 'library_id'])
        .size()
        .reset_index()
    )
    breakpoint_unique_reads.columns = ['cluster_id', 'library_id', 'num_unique_reads']

    breakpoint_library = (
        breakpoint_reads.merge(breakpoint_unique_reads)
        .merge(lib_names)
        .drop(['library_id'], axis=1)
    )

    agg_f = {
        'log_likelihood':np.average,
        'log_cdf':np.average,
        'template_length_1':max,
        'template_length_2':max,
    }

    breakpoint_stats = (
        likelihoods.groupby('cluster_id')
        .agg(agg_f)
        .reset_index()
    )

    breakpoint_stats['template_length_min'] = breakpoint_stats[['template_length_1', 'template_length_2']].min(axis=1)

    breakpoint_counts = (
        likelihoods.groupby('cluster_id')
        .size()
        .reset_index()
    )
    breakpoint_counts.columns = ['cluster_id', 'num_reads']

    breakpoint_unique_counts = (
        likelihoods.drop_duplicates(['cluster_id', 'library_id', 'template_length_1', 'template_length_2'])
        .groupby('cluster_id')
        .size()
        .reset_index()
    )
    breakpoint_unique_counts.columns = ['cluster_id', 'num_unique_reads']

    breakpoints = breakpoints.merge(breakpoint_stats, on='cluster_id', how='inner')
    breakpoints = breakpoints.merge(breakpoint_counts, on='cluster_id', how='inner')
    breakpoints = breakpoints.merge(breakpoint_unique_counts, on='cluster_id', how='inner')

    # Calculate breakpoint type
    def breakpoint_type(row):
        if row['chromosome_1'] != row['chromosome_2']:
            return 'translocation'
        if row['strand_1'] == row['strand_2']:
            return 'inversion'
        positions = sorted([(row['position_{0}'.format(side)], row['strand_{0}'.format(side)]) for side in (1, 2)])
        if positions[0][1] == '+':
            return 'deletion'
        else:
            return 'duplication'

    breakpoints['type'] = breakpoints.apply(breakpoint_type, axis=1)

    # Calculate number inserted at the breakpoint
    def calculate_num_inserted(row):
        if row['inserted'] == '.':
            return 0
        else:
            return len(row['inserted'])

    breakpoints['num_inserted'] = breakpoints.apply(calculate_num_inserted, axis=1)

    # Annotate sequence
    reference_sequences = dict()
    for id, seq in destruct.utils.seq.read_sequences(open(genome_fasta, 'rt')):
        reference_sequences[id] = seq

    breakpoints['sequence'] = breakpoints.apply(lambda row: create_sequence(row, reference_sequences), axis=1)

    # Annotate gene information
    gene_models = pygenes.GeneModels()
    gene_models.load_ensembl_gtf(gtf_filename)

    breakpoints = breakpoints.apply(lambda row: annotate_genes(row, gene_models), axis=1)

    # Annotate database of genomic variants
    dgv = DGVDatabase(dgv_filename)

    breakpoints['dgv_ids'] = breakpoints.apply(lambda row: query_dgv(row, dgv), axis=1)

    breakpoints = breakpoints.rename(columns={'cluster_id':'prediction_id'})

    breakpoints.to_csv(breakpoint_table, sep='\t', na_rep='NA', header=True, index=False)

    breakpoint_library = breakpoint_library.rename(columns={'cluster_id':'prediction_id'})

    breakpoint_library.to_csv(breakpoint_library_table, sep='\t', na_rep='NA', header=True, index=False)



