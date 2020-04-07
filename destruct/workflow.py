import os

import pypeliner
import pypeliner.workflow
import pypeliner.managed as mgd

import destruct.defaultconfig


# Pypeliner contexts
locally = {'local': True}
lowmem = {'mem': 4, 'num_retry': 2, 'mem_retry_factor': 2}
medmem = {'mem': 8, 'num_retry': 2, 'mem_retry_factor': 2}
himem = {'mem': 16, 'num_retry': 2, 'mem_retry_factor': 2}


def create_destruct_workflow(
    bam_filenames,
    breakpoint_table,
    breakpoint_library_table,
    breakpoint_read_table,
    config,
    ref_data_dir,
    raw_data_dir=None,
):
    # Optionally cache raw reads for quicker rerun
    if raw_data_dir is not None:
        mgd_stats = mgd.File(os.path.join(raw_data_dir, '{bylibrary}_stats.txt'), 'bylibrary')
        mgd_reads_1 = mgd.File(os.path.join(raw_data_dir, '{bylibrary}_reads1.fq.gz'), 'bylibrary')
        mgd_reads_2 = mgd.File(os.path.join(raw_data_dir, '{bylibrary}_reads2.fq.gz'), 'bylibrary')
        mgd_sample_1 = mgd.File(os.path.join(raw_data_dir, '{bylibrary}_sample1.fq.gz'), 'bylibrary')
        mgd_sample_2 = mgd.File(os.path.join(raw_data_dir, '{bylibrary}_sample2.fq.gz'), 'bylibrary')

    else:
        mgd_stats = mgd.TempFile('stats.txt', 'bylibrary')
        mgd_reads_1 = mgd.TempFile('reads1.fq.gz', 'bylibrary')
        mgd_reads_2 = mgd.TempFile('reads2.fq.gz', 'bylibrary')
        mgd_sample_1 = mgd.TempFile('sample1.fq.gz', 'bylibrary')
        mgd_sample_2 = mgd.TempFile('sample2.fq.gz', 'bylibrary')

    config = destruct.defaultconfig.get_config(ref_data_dir, config)

    workflow = pypeliner.workflow.Workflow()

    # Set the library ids

    workflow.transform(
        name='create_lib_ids',
        func='destruct.tasks.create_library_ids',
        ret=mgd.TempOutputObj('library_id', 'bylibrary'),
        args=(list(bam_filenames.keys()),)
    )

    # Retrieve discordant reads and stats from bam files

    workflow.commandline(
        name='bamdisc',
        axes=('bylibrary',),
        ctx={'io': 1, 'mem': 8},
        args=(
            'destruct_bamdiscordantfastq',
            '-r',
            '-c', config['bam_max_soft_clipped'],
            '-f', config['bam_max_fragment_length'],
            '-b', mgd.InputFile('bam', 'bylibrary', fnames=bam_filenames),
            '-s', mgd_stats.as_output(),
            '--fastq1', mgd_reads_1.as_output(),
            '--fastq2', mgd_reads_2.as_output(),
            '-t', mgd.TempSpace('bamdisc.tempspace', 'bylibrary'),
            '-n', config['num_read_samples'],
            '--sample1', mgd_sample_1.as_output(),
            '--sample2', mgd_sample_2.as_output(),
        ),
    )

    workflow.subworkflow(
        name='destruct_fastq',
        func='destruct.workflow.create_destruct_fastq_workflow',
        args=(
            mgd_reads_1.as_input(),
            mgd_reads_2.as_input(),
            mgd_sample_1.as_input(),
            mgd_sample_2.as_input(),
            mgd_stats.as_input(),
            mgd.OutputFile(breakpoint_table),
            mgd.OutputFile(breakpoint_library_table),
            mgd.OutputFile(breakpoint_read_table),
            config,
            ref_data_dir,
        ),
        kwargs={
            'raw_data_dir': raw_data_dir,
        },
    )

    return workflow


def create_destruct_fastq_workflow(
    fastq1_filenames,
    fastq2_filenames,
    sample1_filenames,
    sample2_filenames,
    stats_filenames,
    breakpoint_table,
    breakpoint_library_table,
    breakpoint_read_table,
    config,
    ref_data_dir,
    raw_data_dir=None,
):
    workflow = pypeliner.workflow.Workflow()

    # Set the library ids

    workflow.transform(
        name='create_lib_ids',
        func='destruct.tasks.create_library_ids',
        ret=mgd.TempOutputObj('library_id', 'bylibrary'),
        args=(list(fastq1_filenames.keys()),)
    )

    workflow.transform(
        name='readstats',
        axes=('bylibrary',),
        ctx=lowmem,
        func='destruct.tasks.read_stats',
        ret=mgd.TempOutputObj('stats', 'bylibrary'),
        args=(
            mgd.InputFile('stats.txt', 'bylibrary', fnames=stats_filenames),
            config['fragment_length_num_stddevs'],
        ),
    )

    # Align a sample of reads and calculate alignment statistics

    workflow.transform(
        name='prepseed_sample',
        axes=('bylibrary',),
        ctx=medmem,
        func='destruct.tasks.prepare_seed_fastq',
        args=(
            mgd.InputFile('sample1.fq.gz', 'bylibrary', fnames=sample1_filenames),
            mgd.InputFile('sample2.fq.gz', 'bylibrary', fnames=sample2_filenames),
            36,
            mgd.TempOutputFile('sample.seed', 'bylibrary'),
        ),
    )

    workflow.commandline(
        name='bwtrealign_sample',
        axes=('bylibrary',),
        ctx=medmem,
        args=(
            'bowtie',
            config['genome_fasta'],
            mgd.TempInputFile('sample.seed', 'bylibrary'),
            '--chunkmbs', '512',
            '-k', '1000',
            '-m', '1000',
            '--strata',
            '--best',
            '-S',
            '|',
            'destruct_aligntrue',
            '-a', '-',
            '-1', mgd.InputFile('sample1.fq.gz', 'bylibrary', fnames=sample1_filenames),
            '-2', mgd.InputFile('sample2.fq.gz', 'bylibrary', fnames=sample2_filenames),
            '-r', config['genome_fasta'],
            '-g', config['gap_score'],
            '-x', config['mismatch_score'],
            '-m', config['match_score'],
            '--flmin', mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_min'),
            '--flmax', mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_max'),
            '-s', mgd.TempOutputFile('samples.align.true', 'bylibrary'),
        ),
    )

    workflow.transform(
        name='scorestats',
        axes=('bylibrary',),
        ctx=medmem,
        func='destruct.score_stats.create_score_stats',
        args=(
            mgd.TempInputFile('samples.align.true', 'bylibrary'),
            config['match_score'],
            mgd.TempOutputFile('score.stats', 'bylibrary'),
        ),
    )

    # Split discordant fastqs and align

    workflow.transform(
        name='splitfastq1',
        axes=('bylibrary',),
        ctx=lowmem,
        func='destruct.tasks.split_fastq',
        args=(
            mgd.InputFile('reads1.fq.gz', 'bylibrary', fnames=fastq1_filenames),
            int(config['reads_per_split']),
            mgd.TempOutputFile('reads1', 'bylibrary', 'byread'),
        ),
    )

    workflow.transform(
        name='splitfastq2',
        axes=('bylibrary',),
        ctx=lowmem,
        func='destruct.tasks.split_fastq',
        args=(
            mgd.InputFile('reads2.fq.gz', 'bylibrary', fnames=fastq2_filenames),
            int(config['reads_per_split']),
            mgd.TempOutputFile('reads2', 'bylibrary', 'byread', axes_origin=[]),
        ),
    )

    workflow.transform(
        name='prepseed',
        axes=('bylibrary', 'byread'),
        ctx=medmem,
        func='destruct.tasks.prepare_seed_fastq',
        args=(
            mgd.TempInputFile('reads1', 'bylibrary', 'byread'),
            mgd.TempInputFile('reads2', 'bylibrary', 'byread'),
            36,
            mgd.TempOutputFile('reads.seed', 'bylibrary', 'byread'),
        ),
    )

    workflow.commandline(
        name='bwtrealign',
        axes=('bylibrary', 'byread'),
        ctx=medmem,
        args=(
            'bowtie',
            config['genome_fasta'],
            mgd.TempInputFile('reads.seed', 'bylibrary', 'byread'),
            '--chunkmbs', '512',
            '-k', '1000',
            '-m', '1000',
            '--strata',
            '--best',
            '-S',
            '|',
            'destruct_realign2',
            '-l', mgd.TempInputObj('library_id', 'bylibrary'),
            '-a', '-',
            '-1', mgd.TempInputFile('reads1', 'bylibrary', 'byread'),
            '-2', mgd.TempInputFile('reads2', 'bylibrary', 'byread'),
            '-r', config['genome_fasta'],
            '-g', config['gap_score'],
            '-x', config['mismatch_score'],
            '-m', config['match_score'],
            '--flmin', mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_min'),
            '--flmax', mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_max'),
            '--tchimer', config['chimeric_threshold'],
            '--talign', config['alignment_threshold'],
            '--pchimer', config['chimeric_prior'],
            '--tvalid', config['readvalid_threshold'],
            '-z', mgd.TempInputFile('score.stats', 'bylibrary'),
            '--span', mgd.TempOutputFile('spanning.alignments', 'bylibrary', 'byread'),
            '--split', mgd.TempOutputFile('split.alignments', 'bylibrary', 'byread'),
        ),
    )

    workflow.transform(
        name='merge_spanning_1',
        axes=('bylibrary',),
        ctx=lowmem,
        func='destruct.tasks.merge_files_by_line',
        args=(
            mgd.TempInputFile('spanning.alignments', 'bylibrary', 'byread'),
            mgd.TempOutputFile('spanning.alignments_1', 'bylibrary'),
        ),
    )

    workflow.commandline(
        name='filterreads',
        axes=('bylibrary',),
        ctx=lowmem,
        args=(
            'destruct_filterreads',
            '-n', '2',
            '-a', mgd.TempInputFile('spanning.alignments_1', 'bylibrary'),
            '-r', config['satellite_regions'],
            '>', mgd.TempOutputFile('spanning.alignments', 'bylibrary'),
        ),
    )

    workflow.transform(
        name='merge_split_1',
        axes=('bylibrary',),
        ctx=lowmem,
        func='destruct.tasks.merge_files_by_line',
        args=(
            mgd.TempInputFile('split.alignments', 'bylibrary', 'byread'),
            mgd.TempOutputFile('split.alignments', 'bylibrary'),
        ),
    )

    workflow.transform(
        name='merge_spanning_2',
        ctx=lowmem,
        func='destruct.tasks.merge_alignment_files',
        args=(
            mgd.TempInputFile('spanning.alignments', 'bylibrary'),
            mgd.TempOutputFile('spanning.alignments'),
            mgd.TempInputObj('library_id', 'bylibrary'),
        ),
    )

    workflow.transform(
        name='merge_split_2',
        ctx=lowmem,
        func='destruct.tasks.merge_alignment_files',
        args=(
            mgd.TempInputFile('split.alignments', 'bylibrary'),
            mgd.TempOutputFile('split.alignments'),
            mgd.TempInputObj('library_id', 'bylibrary'),
        ),
    )

    # Cluster spanning reads

    workflow.transform(
        name='generate_chrom_args',
        func='destruct.tasks.generate_chromosome_args',
        ret=mgd.TempOutputObj('chrom.args', 'bychromarg'),
        args=(config['chromosomes'],)
    )

    workflow.transform(
        name='write_stats_table',
        ctx=lowmem,
        func='destruct.tasks.write_stats_table',
        args=(
            mgd.TempInputObj('library_id', 'bylibrary'),
            mgd.TempInputObj('stats', 'bylibrary'),
            mgd.TempOutputFile('libstats.tsv'),
        ),
    )

    workflow.commandline(
        name='cluster',
        axes=('bychromarg',),
        ctx=medmem,
        args=(
            'destruct_mclustermatepairs',
            '-a', mgd.TempInputFile('spanning.alignments'),
            '-s', mgd.TempInputFile('libstats.tsv'),
            '-c', mgd.TempOutputFile('clusters', 'bychromarg'),
            mgd.TempInputObj('chrom.args', 'bychromarg'),
            '--clustmin', config['cluster_readcount_threshold'],
            '--fragmax', config['fragment_length_max'],
        ),
    )
    
    # Predict breakpoints from split reads

    workflow.transform(
        name='predict_breaks',
        axes=('bychromarg',),
        ctx=medmem,
        func='destruct.predict_breaks.predict_breaks',
        args=(
            mgd.TempInputFile('clusters', 'bychromarg'),
            mgd.TempInputFile('spanning.alignments'),
            mgd.TempInputFile('split.alignments'),
            mgd.TempOutputFile('breakpoints_2', 'bychromarg'),
        ),
    )

    workflow.transform(
        name='merge_clusters',
        ctx=lowmem,
        func='destruct.tasks.merge_clusters',
        args=(
            mgd.TempInputFile('clusters', 'bychromarg'),
            mgd.TempInputFile('breakpoints_2', 'bychromarg'),
            mgd.TempOutputFile('clusters'),
            mgd.TempOutputFile('breakpoints_2'),
            mgd.TempOutputFile('merge_clusters.debug'),
        ),
    )

    # Realign reads to breakpoints

    workflow.commandline(
        name='realigntobreaks',
        axes=('bylibrary', 'byread'),
        ctx=medmem,
        args=(
            'destruct_realigntobreaks2',
            '-r', config['genome_fasta'],
            '-b', mgd.TempInputFile('breakpoints_2'),
            '-c', mgd.TempInputFile('clusters'),
            '-g', config['gap_score'],
            '-x', config['mismatch_score'],
            '-m', config['match_score'],
            '--flmax', mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_max'),
            '--span', mgd.TempInputFile('spanning.alignments', 'bylibrary', 'byread'),
            '-1', mgd.TempInputFile('reads1', 'bylibrary', 'byread'),
            '-2', mgd.TempInputFile('reads2', 'bylibrary', 'byread'),
            '--realignments', mgd.TempOutputFile('realignments', 'bylibrary', 'byread'),
        ),
    )

    # Calculate likelihoods based on realignments

    workflow.transform(
        name='calculate_realignment_likelihoods',
        axes=('bylibrary', 'byread'),
        ctx=medmem,
        func='destruct.predict_breaks.calculate_realignment_likelihoods',
        args=(
            mgd.TempInputFile('breakpoints_2'),
            mgd.TempInputFile('realignments', 'bylibrary', 'byread'),
            mgd.TempInputFile('score.stats', 'bylibrary'),
            mgd.TempOutputFile('likelihoods_2', 'bylibrary', 'byread'),
            config['match_score'],
            mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_mean'),
            mgd.TempInputObj('stats', 'bylibrary').prop('fragment_length_stddev'),
        ),
    )

    workflow.transform(
        name='merge_likelihoods_1',
        axes=('bylibrary',),
        ctx=lowmem,
        func='destruct.tasks.merge_sorted_files_by_line',
        args=(
            mgd.TempInputFile('likelihoods_2', 'bylibrary', 'byread'),
            mgd.TempOutputFile('likelihoods_2', 'bylibrary'),
            mgd.TempSpace('merge_likelihoods_1_temp', 'bylibrary'),
            '1',
        ),
    )

    workflow.transform(
        name='merge_likelihoods_2',
        ctx=lowmem,
        func='destruct.tasks.merge_sorted_files_by_line',
        args=(
            mgd.TempInputFile('likelihoods_2', 'bylibrary'),
            mgd.TempOutputFile('likelihoods_2'),
            mgd.TempSpace('merge_likelihoods_2_temp'),
            '1',
        ),
    )

    # Set cover for multi mapping reads

    workflow.transform(
        name='calc_weights',
        ctx=medmem,
        func='destruct.predict_breaks.calculate_cluster_weights',
        args=(
            mgd.TempInputFile('breakpoints_2'),
            mgd.TempOutputFile('cluster_weights'),
        ),
    )

    workflow.commandline(
        name='setcover',
        ctx=medmem,
        args=(
            'destruct_setcover',
            '-c', mgd.TempInputFile('clusters'),
            '-w', mgd.TempInputFile('cluster_weights'),
            '-a', mgd.TempOutputFile('clusters_setcover'),
        ),
    )

    # Select cluster based on setcover

    workflow.transform(
        name='select_clusters',
        ctx=medmem,
        func='destruct.predict_breaks.select_clusters',
        args=(
            mgd.TempInputFile('clusters_setcover'),
            mgd.TempInputFile('breakpoints_2'),
            mgd.TempOutputFile('breakpoints_1'),
            mgd.TempInputFile('likelihoods_2'),
            mgd.TempOutputFile('likelihoods_1'),
        ),
    )

    # Select prediction based on max likelihood

    workflow.transform(
        name='select_predictions',
        ctx=himem,
        func='destruct.predict_breaks.select_predictions',
        args=(
            mgd.TempInputFile('breakpoints_1'),
            mgd.TempOutputFile('breakpoints'),
            mgd.TempInputFile('likelihoods_1'),
            mgd.TempOutputFile('likelihoods'),
            config['mate_score_threshold'],
            config['template_length_min_threshold'],
            config['min_alignment_log_likelihood'],
        ),
    )

    # Optionally tabulate supporting reads

    workflow.transform(
        name='tabreads',
        ctx=medmem,
        func='destruct.tasks.tabulate_reads',
        args=(
            mgd.TempInputFile('clusters_setcover'),
            mgd.TempInputFile('likelihoods'),
            mgd.TempInputObj('library_id', 'bylibrary'),
            mgd.InputFile('reads1.fq.gz', 'bylibrary', fnames=fastq1_filenames),
            mgd.InputFile('reads2.fq.gz', 'bylibrary', fnames=fastq2_filenames),
            mgd.TempOutputFile('breakreads.table.unsorted'),
        ),
    )

    workflow.commandline(
        name='sortreads',
        ctx=medmem,
        args=(
            'sort', '-n',
            mgd.TempInputFile('breakreads.table.unsorted'),
            '>', mgd.OutputFile(breakpoint_read_table),
        ),
    )


    # Tabulate results

    workflow.transform(
        name='tabulate',
        ctx=himem,
        func='destruct.tasks.tabulate_results',
        args=(
            mgd.TempInputFile('breakpoints'),
            mgd.TempInputFile('likelihoods'),
            mgd.TempInputObj('library_id', 'bylibrary'),
            config['genome_fasta'],
            config['gtf_filename'],
            config['dgv_filename'],
            mgd.OutputFile(breakpoint_table),
            mgd.OutputFile(breakpoint_library_table),
        ),
    )

    return workflow

