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
        name='aln1',
        axes=('byread',),
        ctx=himem,
        args=(
            'bwa', 'aln',
            genome_fasta_filename,
            mgd.TempInputFile('fastq1', 'byread'),
            '>',
            mgd.TempOutputFile('sai1', 'byread'),
        ),
    )

    workflow.commandline(
        name='aln2',
        axes=('byread',),
        ctx=himem,
        args=(
            'bwa', 'aln',
            genome_fasta_filename,
            mgd.TempInputFile('fastq2', 'byread'),
            '>',
            mgd.TempOutputFile('sai2', 'byread'),
        ),
    )

    workflow.commandline(
        name='sampe',
        axes=('byread',),
        ctx=himem,
        args=(
            'bwa', 'sampe',
            '-r', read_group_str,
            genome_fasta_filename,
            mgd.TempInputFile('sai1', 'byread'),
            mgd.TempInputFile('sai2', 'byread'),
            mgd.TempInputFile('fastq1', 'byread'),
            mgd.TempInputFile('fastq2', 'byread'),
            '>',
            mgd.TempOutputFile('sam', 'byread'),
        ),
    )

    workflow.commandline(
        name='create_bam',
        axes=('byread',),
        ctx=lowmem,
        args=(
            'samtools', 'view', '-bt',
            genome_fasta_filename+'.fai',
            mgd.TempInputFile('sam', 'byread'),
            '>',
            mgd.TempOutputFile('bam', 'byread'),
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


