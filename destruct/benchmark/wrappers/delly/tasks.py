import os
import pysam
import pandas as pd

import pypeliner.commandline


def _rename_index(out_file):
    if out_file.endswith('.tmp'):
        out_index = out_file + '.csi'
        renamed_out_index = out_file[:-4] + '.csi'
        try:
            os.remove(renamed_out_index)
        except OSError:
            pass
        os.rename(out_index, renamed_out_index)


def _write_empty_vcf(out_file, bam_files):
    with open(out_file+'.vcf', 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')
        for bam_file in bam_files:
            sample_id = os.path.splitext(os.path.basename(bam_file))[0]
            f.write('\t{}'.format(sample_id))
        f.write('\n')

    pypeliner.commandline.execute('bcftools', 'convert', '-O', 'b', out_file+'.vcf', '-o', out_file)
    pypeliner.commandline.execute('bcftools', 'index', out_file)
    os.remove(out_file+'.vcf')


def get_sv_types(ref_genome_fasta_file):

    # Check for a genome with a single chromosome
    num_chromosomes = 0
    with open(ref_genome_fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                num_chromosomes += 1
            if num_chromosomes > 1:
                break

    # Translocations only for genomes with multiple chromosomes
    if num_chromosomes > 1:
        return ('DEL', 'DUP', 'INV', 'TRA', 'INS')
    else:
        return ('DEL', 'DUP', 'INV', 'INS')


def run_delly_call(sv_type, delly_excl_chrom, ref_genome_fasta_file, bam_files, out_file):
    delly_args = [
        'delly', 'call',
        '-t', sv_type,
        '-x', delly_excl_chrom,
        '-g', ref_genome_fasta_file,
        '-o', out_file,
    ]

    delly_args += bam_files
    
    pypeliner.commandline.execute(*delly_args)

    if not os.path.exists(out_file):
        _write_empty_vcf(out_file, bam_files)

    _rename_index(out_file)


def run_delly_filter(sv_type, sample_ids, control_id, sample_file, ref_genome_fasta_file, in_file, out_file):
    with open(sample_file, 'w') as f:
        for sample_id in sample_ids:
            if control_id == sample_id:
                sample_type = 'control'
            else:
                sample_type = 'tumor'
            f.write('{}\t{}\n'.format(sample_id, sample_type))

    delly_args = [
        'delly', 'filter',
        '-t', sv_type,
        '-f', 'somatic',
        '-o', out_file,
        '-s', sample_file,
        '-g', ref_genome_fasta_file,
        in_file,
    ]

    pypeliner.commandline.execute(*delly_args)
    _rename_index(out_file)


def convert_vcf(bcf_filename, output_filename, control_id=None):
    bcf_reader = pysam.VariantFile(bcf_filename, 'rb')

    breakpoint_table = list()
    breakpoint_library_table = list()

    for row in bcf_reader:
        prediction_id = row.id

        chrom_1 = row.chrom
        chrom_2 = row.info['CHR2']

        strand_1, strand_2 = [('-', '+')[a == '3'] for a in row.info['CT'].split('to')]

        coord_1 = row.pos
        coord_2 = row.info['END']

        if 'LowQual' in row.filter:
            qual = 0
        else:
            qual = 1

        breakpoint_table.append((prediction_id, chrom_1, chrom_2, strand_1, strand_2, coord_1, coord_2, qual))

        for sample, call in row.samples.items():
            num_spanning = call['DV']
            num_split = call['RV']

            breakpoint_library_table.append((prediction_id, sample, num_spanning, num_split))

    breakpoint_table = pd.DataFrame(
        breakpoint_table,
        columns=[
            'prediction_id',
            'chromosome_1', 'chromosome_2',
            'strand_1', 'strand_2',
            'position_1', 'position_2',
            'qual',
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

    if len(breakpoint_table.index) == 0:
        raise Exception('No variants')

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


