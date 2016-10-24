import os

import pypeliner.workflow
import pypeliner.managed as mgd

import destruct.benchmark.align.bwa.tasks


def bwa_align_workflow(
    genome_fasta_filename,
    fastq_1_filename,
    fastq_2_filename,
    bam_filename,
    reads_per_job=int(1e6),
    read_group_str=None,
):
    if read_group_str is None:
        read_group_str = '@RG\\tID:S1\\tSM:sample_1'

    lowmem = {'mem':1}
    himem = {'mem':8}

    if not os.path.exists(genome_fasta_filename+'.bwt'):
        raise Exception('No index for ' + genome_fasta_filename)

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='split1',
        ctx=lowmem,
        func=destruct.benchmark.align.bwa.tasks.split_fastq,
        args=(
            mgd.InputFile(fastq_1_filename),
            reads_per_job,
            mgd.TempOutputFile('fastq1', 'byread'),
        ),
    )

    workflow.transform(
        name='split2',
        ctx=lowmem,
        func=destruct.benchmark.align.bwa.tasks.split_fastq,
        args=(
            mgd.InputFile(fastq_2_filename), 
            reads_per_job, 
            mgd.TempOutputFile('fastq2', 'byread', axes_origin=[]),
        ),
    )

    workflow.commandline(
        name='align',
        axes=('byread',),
        ctx=himem,
        args=(
            'bwa', 'mem',
            '-R', read_group_str,
            genome_fasta_filename,
            mgd.TempInputFile('fastq1', 'byread'),
            mgd.TempInputFile('fastq2', 'byread'),
            '|',
            'samtools', 'view', '-bt',
            genome_fasta_filename+'.fai',
            '-',
            '-o', mgd.TempOutputFile('bam', 'byread'),
        ),
    )

    workflow.transform(
        name='merge_bams',
        ctx=lowmem,
        func=destruct.benchmark.align.bwa.tasks.merge_bam,
        args=(
            mgd.TempInputFile('bam', 'byread'),
            mgd.OutputFile(bam_filename),
        ),
    )

    return workflow


