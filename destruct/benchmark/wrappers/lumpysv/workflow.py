import pypeliner
import pypeliner.managed as mgd

import destruct.workflow
import destruct.benchmark.wrappers.utils
import destruct.benchmark.wrappers.lumpysv.tasks


def create_lumpysv_wrapper_workflow(bam_filenames, output_filename, raw_data_dir, control_id=None):
    bams = list()
    for lib_id, bam_filename in bam_filenames.items():
        bams += [destruct.benchmark.wrappers.utils.symlink(bam_filename, link_name='{0}.bam'.format(lib_id), link_directory=raw_data_dir)]
        destruct.benchmark.wrappers.utils.symlink(bam_filename+'.bai', link_name='{0}.bam.bai'.format(lib_id), link_directory=raw_data_dir)

    workflow = pypeliner.workflow.Workflow(default_ctx={'mem': 4})

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=list(bam_filenames.keys()),
    )

    workflow.commandline(
        name='extract_discordants',
        axes=('sample_id',),
        args=(
            'samtools', 'view', '-b', '-F', '1294',
            mgd.InputFile('bam', 'sample_id', fnames=bam_filenames),
            '>',
            mgd.TempOutputFile('discordants.unsorted.bam', 'sample_id'),
        ),
    )

    workflow.commandline(
        name='extract_splitters',
        axes=('sample_id',),
        args=(
            'samtools', 'view', '-h',
            mgd.InputFile('bam', 'sample_id', fnames=bam_filenames),
            '|',
            'lumpy_extractSplitReads_BwaMem', '-i', 'stdin',
            '|',
            'samtools', 'view', '-Sb', '-',
            '>',
            mgd.TempOutputFile('splitters.unsorted.bam', 'sample_id'),
        ),
    )

    workflow.commandline(
        name='sort_discordants',
        axes=('sample_id',),
        args=(
            'samtools', 'sort',
            mgd.TempInputFile('discordants.unsorted.bam', 'sample_id'),
            '-o',
            mgd.TempOutputFile('discordants.bam', 'sample_id'),
        ),
    )

    workflow.commandline(
        name='sort_splitters',
        axes=('sample_id',),
        args=(
            'samtools', 'sort',
            mgd.TempInputFile('splitters.unsorted.bam', 'sample_id'),
            '-o',
            mgd.TempOutputFile('splitters.bam', 'sample_id'),
        ),
    )

    workflow.transform(
        name='run_lumpyexpress',
        func=destruct.benchmark.wrappers.lumpysv.tasks.run_lumpyexpress,
        args=(
            mgd.InputFile('bam', 'sample_id', fnames=bam_filenames),
            mgd.TempInputFile('splitters.bam', 'sample_id'),
            mgd.TempInputFile('discordants.bam', 'sample_id'),
            mgd.TempOutputFile('results.vcf'),
        ),
    )

    workflow.transform(
        name='vcf_to_bcf',
        func=destruct.benchmark.wrappers.lumpysv.tasks.vcf_to_bcf,
        args=(
            mgd.TempInputFile('results.vcf'),
            mgd.TempOutputFile('results.bcf'),
        ),
    )

    workflow.transform(
        name='convert_bcf',
        func=destruct.benchmark.wrappers.lumpysv.tasks.convert_bcf,
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
        args=(
            mgd.TempInputFile('results.bcf'),
            mgd.OutputFile(output_filename),
        ),
        kwargs={
            'control_id': control_id,
        },
    )

    return workflow


def setup_lumpysv(test_config, **kwargs):
    pass
