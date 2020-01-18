import os
import pysam
import pandas as pd

import pypeliner.commandline


def run_lumpyexpress(bam_filenames, splitters_filenames, discordants_filenames, results_vcf):
    bam_arg = []
    splitters_arg = []
    discordants_arg = []
    for sample_id in bam_filenames.keys():
        bam_arg.append(bam_filenames[sample_id])
        splitters_arg.append(splitters_filenames[sample_id])
        discordants_arg.append(discordants_filenames[sample_id])

    bam_arg = ','.join(bam_arg)
    splitters_arg = ','.join(splitters_arg)
    discordants_arg = ','.join(discordants_arg)

    pypeliner.commandline.execute(
        'lumpyexpress',
        '-B', bam_arg,
        '-S', splitters_arg,
        '-D', discordants_arg,
        '-o', results_vcf,
    )


def vcf_to_bcf(vcf_filename, bcf_filename):
    sorted_vcf = vcf_filename + '.sorted'
    pypeliner.commandline.execute('vcf-sort', vcf_filename, '>', sorted_vcf)

    gz_vcf = vcf_filename + '.sorted.gz'
    pypeliner.commandline.execute('bgzip', sorted_vcf)

    gz_vcf_tbi = vcf_filename + '.sorted.gz.tbi'
    pypeliner.commandline.execute('tabix', gz_vcf)

    pypeliner.commandline.execute('bcftools', 'convert', '-O', 'b', gz_vcf, '-o', bcf_filename)

    os.remove(gz_vcf)
    os.remove(gz_vcf_tbi)


def convert_bcf(bcf_filename, output_filename, control_id=None):
    bcf_reader = pysam.VariantFile(bcf_filename, 'rb')

    breakpoint_table = list()
    breakpoint_library_table = list()

    for row in bcf_reader:
        if row.info.get('SECONDARY', False):
            continue

        prediction_id = row.id

        assert row.info['SVTYPE'] in ('DEL', 'DUP', 'INV', 'BND')

        chrom_1 = row.chrom
        coord_1 = row.pos
        if row.info['SVTYPE'] == 'BND':
            if '[' in row.alts[0]:
                chrom_2, coord_2 = row.alts[0].split('[')[1].split(':')
                coord_2 = int(coord_2)
            elif ']' in row.alts[0]:
                chrom_2, coord_2 = row.alts[0].split(']')[1].split(':')
                coord_2 = int(coord_2)
        else:
            chrom_2 = row.chrom
            coord_2 = row.stop
            assert coord_1 < coord_2

        breakpoint_ids = []
        assert isinstance(row.info['STRANDS'], tuple)
        for strand_info, breakpoint_id in zip(row.info['STRANDS'], 'AB'):
            (strand_1, strand_2), breakpoint_read_count = strand_info.split(':')

            if row.info['SVTYPE'] == 'DEL':
                assert (strand_1, strand_2) == ('+', '-')

            elif row.info['SVTYPE'] == 'DUP':
                assert (strand_1, strand_2) == ('-', '+')

            elif row.info['SVTYPE'] == 'DUP':
                assert strand_1 == strand_2

            svtype = row.info['SVTYPE']

            breakpoint_table.append((prediction_id + '_' + breakpoint_id, chrom_1, chrom_2, strand_1, strand_2, coord_1, coord_2, svtype, breakpoint_read_count))
            breakpoint_ids.append(prediction_id + '_' + breakpoint_id)

        for sample, call in row.samples.items():
            num_spanning = call['PE']
            num_split = call['SR']

            for breakpoint_id in breakpoint_ids:
                breakpoint_library_table.append((breakpoint_id, sample, num_spanning, num_split))

    breakpoint_table = pd.DataFrame(
        breakpoint_table,
        columns=[
            'prediction_id',
            'chromosome_1', 'chromosome_2',
            'strand_1', 'strand_2',
            'position_1', 'position_2',
            'svtype',
            'breakpoint_read_count',
        ]
    )

    breakpoint_library_table = pd.DataFrame(
        breakpoint_library_table,
        columns=[
            'prediction_id',
            'library',
            'num_spanning',
            'num_split'
        ]
    )

    spanning_read_counts = breakpoint_library_table.groupby('prediction_id')['num_spanning'].sum().reset_index()
    split_read_counts = breakpoint_library_table.groupby('prediction_id')['num_split'].sum().reset_index()

    total_read_counts = breakpoint_library_table.set_index(['prediction_id', 'library']).sum(axis=1).unstack()
    total_read_counts.columns = [a + '_count' for a in total_read_counts.columns]
    total_read_counts = total_read_counts.reset_index()

    breakpoint_table = breakpoint_table.merge(spanning_read_counts, on='prediction_id')
    breakpoint_table = breakpoint_table.merge(split_read_counts, on='prediction_id')
    breakpoint_table = breakpoint_table.merge(total_read_counts, on='prediction_id')

    # Filter based on evidence in control dataset
    if control_id is not None:
         breakpoint_table = breakpoint_table[breakpoint_table['{0}_count'.format(control_id)] == 0]          

    breakpoint_table.to_csv(output_filename, sep='\t', index=False)


