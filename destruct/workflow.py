import os

import pypeliner
import pypeliner.workflow
import pypeliner.managed as mgd

import destruct.score_stats
import destruct.utils.misc
import destruct.utils.seq
import destruct.predict_breaks
import destruct.tasks
import destruct.defaultconfig


# Pypeliner contexts
locally = {'local': True}
lowmem = {'mem': 4, 'num_retry': 2, 'mem_retry_factor': 2}
medmem = {'mem': 8, 'num_retry': 2, 'mem_retry_factor': 2}
himem = {'mem': 16, 'num_retry': 2, 'mem_retry_factor': 2}


def align(
    stats_filename,
    reads1_filename,
    reads2_filename,
    sample1_filename,
    sample2_filename,
    spanning_filename,
    split_filename,
    score_stats_filename,
    library_id,
    config,
):
    workflow = pypeliner.workflow.Workflow()

    # Retrieve discordant reads and stats from bam files
    workflow.transform(
        name='readstats',
        ctx=lowmem,
        func=destruct.tasks.read_stats,
        ret=mgd.TempOutputObj('stats'),
        args=(
            mgd.InputFile(stats_filename),
            config['fragment_length_num_stddevs'],
        ),
    )

    # Align a sample of reads and calculate alignment statistics
    workflow.transform(
        name='prepseed_sample',
        ctx=medmem,
        func=destruct.tasks.prepare_seed_fastq,
        args=(
            mgd.InputFile(sample1_filename),
            mgd.InputFile(sample2_filename),
            36,
            mgd.TempOutputFile('sample.seed'),
        ),
    )

    workflow.commandline(
        name='bwtrealign_sample',
        ctx=medmem,
        args=(
            'bowtie',
            config['genome_fasta'],
            mgd.TempInputFile('sample.seed'),
            '--chunkmbs', '512',
            '-k', '1000',
            '-m', '1000',
            '--strata',
            '--best',
            '-S',
            '|',
            'destruct_aligntrue',
            '-a', '-',
            '-1', mgd.InputFile(sample1_filename),
            '-2', mgd.InputFile(sample2_filename),
            '-r', config['genome_fasta'],
            '-g', config['gap_score'],
            '-x', config['mismatch_score'],
            '-m', config['match_score'],
            '--flmin', mgd.TempInputObj('stats').prop('fragment_length_min'),
            '--flmax', mgd.TempInputObj('stats').prop('fragment_length_max'),
            '-s', mgd.TempOutputFile('samples.align.true'),
        ),
    )

    workflow.transform(
        name='scorestats',
        ctx=medmem,
        func=destruct.score_stats.create_score_stats,
        args=(
            mgd.TempInputFile('samples.align.true'),
            config['match_score'],
            mgd.OutputFile(score_stats_filename),
        ),
    )

    # Split discordant fastqs and align
    workflow.transform(
        name='splitfastq1',
        ctx=lowmem,
        func=destruct.tasks.split_fastq,
        args=(
            mgd.InputFile(reads1_filename),
            int(config['reads_per_split']),
            mgd.TempOutputFile('reads1', 'byread'),
        ),
    )

    workflow.transform(
        name='splitfastq2',
        ctx=lowmem,
        func=destruct.tasks.split_fastq,
        args=(
            mgd.InputFile(reads2_filename),
            int(config['reads_per_split']),
            mgd.TempOutputFile('reads2', 'byread', axes_origin=[]),
        ),
    )

    workflow.transform(
        name='prepseed',
        axes=('byread',),
        ctx=medmem,
        func=destruct.tasks.prepare_seed_fastq,
        args=(
            mgd.TempInputFile('reads1', 'byread'),
            mgd.TempInputFile('reads2', 'byread'),
            36,
            mgd.TempOutputFile('reads.seed', 'byread'),
        ),
    )

    workflow.commandline(
        name='bwtrealign',
        axes=('byread',),
        ctx=medmem,
        args=(
            'bowtie',
            config['genome_fasta'],
            mgd.TempInputFile('reads.seed', 'byread'),
            '--chunkmbs', '512',
            '-k', '1000',
            '-m', '1000',
            '--strata',
            '--best',
            '-S',
            '|',
            'destruct_realign2',
            '-l', library_id,
            '-a', '-',
            '-1', mgd.TempInputFile('reads1', 'byread'),
            '-2', mgd.TempInputFile('reads2', 'byread'),
            '-r', config['genome_fasta'],
            '-g', config['gap_score'],
            '-x', config['mismatch_score'],
            '-m', config['match_score'],
            '--flmin', mgd.TempInputObj('stats').prop('fragment_length_min'),
            '--flmax', mgd.TempInputObj('stats').prop('fragment_length_max'),
            '--tchimer', config['chimeric_threshold'],
            '--talign', config['alignment_threshold'],
            '--pchimer', config['chimeric_prior'],
            '--tvalid', config['readvalid_threshold'],
            '-z', mgd.InputFile(score_stats_filename),
            '--span', mgd.TempOutputFile('spanning.alignments', 'byread'),
            '--split', mgd.TempOutputFile('split.alignments', 'byread'),
        ),
    )

    workflow.transform(
        name='merge_spanning_1',
        ctx=lowmem,
        func=destruct.tasks.merge_files_by_line,
        args=(
            mgd.TempInputFile('spanning.alignments', 'byread'),
            mgd.OutputFile(spanning_filename),
        ),
    )

    workflow.transform(
        name='merge_split_1',
        axes=('bylibrary',),
        ctx=lowmem,
        func=destruct.tasks.merge_files_by_line,
        args=(
            mgd.TempInputFile('split.alignments', 'byread'),
            mgd.OutputFile(split_filename),
        ),
    )

    return workflow


def cluster(
    spanning_alignments_filename,
    split_alignments_filename,
    clusters_filename,
    breakpoints_filename,
):
    workflow = pypeliner.workflow.Workflow()

    # Cluster spanning reads
    workflow.setobj(
        obj=mgd.TempOutputObj('chrom.args', 'bychromarg'),
        value=destruct.tasks.generate_chromosome_args(config['chromosomes']),
    )

    workflow.transform(
        name='write_stats_table',
        ctx=lowmem,
        func=destruct.tasks.write_stats_table,
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
            '-a', mgd.InputFile(spanning_alignments_filename),
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
        func=destruct.predict_breaks.predict_breaks,
        args=(
            mgd.TempInputFile('clusters', 'bychromarg'),
            mgd.TempInputFile(spanning_alignments_filename),
            mgd.TempInputFile(split_alignments_filename),
            mgd.TempOutputFile('breakpoints_2', 'bychromarg'),
        ),
    )

    workflow.transform(
        name='merge_clusters',
        ctx=lowmem,
        func=destruct.tasks.merge_clusters,
        args=(
            mgd.TempInputFile('clusters', 'bychromarg'),
            mgd.TempInputFile('breakpoints_2', 'bychromarg'),
            mgd.TempOutputFile(clusters_filename),
            mgd.TempOutputFile(breakpoints_filename),
            mgd.TempOutputFile('merge_clusters.debug'),
        ),
    )

    return workflow


def realign(
    reads1_filename,
    reads2_filename,
    clusters_filename,
    breakpoints_filename,
    likelihoods_filename,
    spanning_alignments_filename,
    score_stats_filename,
    stats,
    library_id,
    config,
):
    workflow = pypeliner.workflow.Workflow()

    # Split discordant fastqs and align
    workflow.transform(
        name='splitfastq1',
        ctx=lowmem,
        func=destruct.tasks.split_fastq,
        args=(
            mgd.InputFile(reads1_filename),
            int(config['reads_per_split']),
            mgd.TempOutputFile('reads1', 'byread'),
        ),
    )

    workflow.transform(
        name='splitfastq2',
        ctx=lowmem,
        func=destruct.tasks.split_fastq,
        args=(
            mgd.InputFile(reads2_filename),
            int(config['reads_per_split']),
            mgd.TempOutputFile('reads2', 'byread', axes_origin=[]),
        ),
    )

    # Realign reads to breakpoints
    workflow.commandline(
        name='realigntobreaks',
        axes=('byread',),
        ctx=medmem,
        args=(
            'destruct_realigntobreaks2',
            '-r', config['genome_fasta'],
            '-b', mgd.InputFile(breakpoints_filename),
            '-c', mgd.InputFile(clusters_filename),
            '-g', config['gap_score'],
            '-x', config['mismatch_score'],
            '-m', config['match_score'],
            '--flmax', stats.fragment_length_max,
            '--span', mgd.InputFile(spanning_alignments_filename),
            '-1', mgd.TempInputFile('reads1', 'byread'),
            '-2', mgd.TempInputFile('reads2', 'byread'),
            '--realignments', mgd.TempOutputFile('realignments', 'byread', axes_origin=[]),
            '-l', library_id,
        ),
    )

    # Calculate likelihoods based on realignments
    workflow.transform(
        name='calculate_realignment_likelihoods',
        axes=('byread',),
        ctx=medmem,
        func=destruct.predict_breaks.calculate_realignment_likelihoods,
        args=(
            mgd.InputFile(breakpoints_filename),
            mgd.TempInputFile('realignments', 'byread'),
            mgd.InputFile(score_stats_filename),
            mgd.TempOutputFile('likelihoods_2', 'byread'),
            config['match_score'],
            stats.fragment_length_mean,
            stats.fragment_length_stddev,
        ),
    )

    workflow.transform(
        name='merge_likelihoods_1',
        axes=('',),
        ctx=lowmem,
        func=destruct.tasks.merge_sorted_files_by_line,
        args=(
            mgd.TempInputFile('likelihoods_2', 'byread'),
            mgd.OutputFile(likelihoods_filename),
            mgd.TempSpace('merge_likelihoods_1_temp'),
            '1',
        ),
    )

    return workflow


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
    workflow.setobj(
        obj=mgd.TempOutputObj('library_id', 'bylibrary'),
        value=destruct.tasks.create_library_ids(bam_filenames.keys()),
    )

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

    # Retrieve discordant reads and stats from bam files
    workflow.transform(
        name='readstats',
        axes=('bylibrary',),
        ctx=lowmem,
        func=destruct.tasks.read_stats,
        ret=mgd.TempOutputObj('stats', 'bylibrary'),
        args=(
            mgd_stats.as_input(),
            config['fragment_length_num_stddevs'],
        ),
    )

    workflow.subworkflow(
        name='align',
        axes=('bylibrary',),
        func=align,
        args=(
            mgd_stats.as_input(),
            mgd_reads_1.as_input(),
            mgd_reads_2.as_input(),
            mgd_sample_1.as_input(),
            mgd_sample_2.as_input(),
            mgd.TempOutputFile('spanning.alignments_1', 'bylibrary'),
            mgd.TempOutputFile('split.alignments', 'bylibrary'),
            mgd.TempOutputFile('score.stats', 'bylibrary'),
            mgd.TempInputObj('library_id', 'bylibrary'),
            config,
        )
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
        name='merge_spanning_2',
        ctx=lowmem,
        func=destruct.tasks.merge_alignment_files,
        args=(
            mgd.TempInputFile('spanning.alignments', 'bylibrary'),
            mgd.TempOutputFile('spanning.alignments'),
            mgd.TempInputObj('library_id', 'bylibrary'),
        ),
    )

    workflow.transform(
        name='merge_split_2',
        ctx=lowmem,
        func=destruct.tasks.merge_alignment_files,
        args=(
            mgd.TempInputFile('split.alignments', 'bylibrary'),
            mgd.TempOutputFile('split.alignments'),
            mgd.TempInputObj('library_id', 'bylibrary'),
        ),
    )

    workflow.subworkflow(

    spanning_alignments_filename,
    split_alignments_filename,
    clusters_filename,
    breakpoints_filename,
    
    )

    workflow.subworkflow(
        name='realign',
        func=realign,
        axes=('bylibrary',),
        args=(
            mgd_reads_1.as_input(),
            mgd_reads_2.as_input(),
            mgd.TempInputFile('clusters'),
            mgd.TempInputFile('breakpoints_2'),
            mgd.TempOutputFile('likelihoods_2', 'bylibrary'),
            mgd.TempInputFile('spanning.alignments', 'bylibrary'),
            mgd.TempInputFile('score.stats', 'bylibrary'),
            mgd.TempInputObj('stats', 'bylibrary'),
            mgd.TempInputObj('library_id', 'bylibrary'),
            config,
        ),
    )

    workflow.transform(
        name='merge_likelihoods_2',
        ctx=lowmem,
        func=destruct.tasks.merge_sorted_files_by_line,
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
        func=destruct.predict_breaks.calculate_cluster_weights,
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
        func=destruct.predict_breaks.select_clusters,
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
        func=destruct.predict_breaks.select_predictions,
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
        func=destruct.tasks.tabulate_reads,
        args=(
            mgd.TempInputFile('clusters_setcover'),
            mgd.TempInputObj('library_id', 'bylibrary'),
            mgd_reads_1.as_input(),
            mgd_reads_2.as_input(),
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
        func=destruct.tasks.tabulate_results,
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
