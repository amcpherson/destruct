import os
import argparse
import yaml

import pypeliner
import pypeliner.workflow
import pypeliner.managed as mgd

import destruct.benchmark.align.bwa.workflow
import destruct.benchmark.destruct_test
import destruct.benchmark.create_breakpoint_simulation


if __name__ == '__main__':

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('simconfig',
                           help='Simulation configuration filename')

    argparser.add_argument('bam',
                           help='Source bam filename')

    argparser.add_argument('ref',
                           help='Reference genome for source bam')

    argparser.add_argument('outdir',
                           help='Output directory')

    argparser.add_argument('tool_defs',
                           help='Tool Definition Filename')

    argparser.add_argument('--config',
                           help='Configuration filename')

    argparser.add_argument('--tool_names', nargs='+',
                           help='Tools to benchmark')

    args = vars(argparser.parse_args())

    config = {}

    if args['config'] is not None:
        execfile(args['config'], {}, config)

    config.update(args)

    pyp = pypeliner.app.Pypeline(config=config)

    tool_defs = yaml.load(open(args['tool_defs']))
    del tool_defs['databases']

    try:
        os.makedirs(args['outdir'])
    except OSError:
        pass

    ctx = {'mem':4}

    workflow = pypeliner.workflow.Workflow(default_ctx=ctx)

    workflow.transform(
        name='read_params',
        func=destruct.benchmark.destruct_test.read_simulation_params,
        ret=mgd.TempOutputObj('simulation.params'),
        args=(mgd.InputFile(args['simconfig']),),
    )

    workflow.transform(
        name='create_sim',
        func=destruct.benchmark.create_breakpoint_simulation.create_breakpoints,
        args=(
            mgd.TempInputObj('simulation.params'),
            mgd.InputFile(args['ref']),
            mgd.OutputFile(os.path.join(args['outdir'], 'simulated.fasta')),
            mgd.OutputFile(os.path.join(args['outdir'], 'simulated.tsv')),
        ),
    )

    workflow.transform(
        name='partition',
        func=destruct.benchmark.destruct_test.partition_bam,
        args=(
            mgd.InputFile(args['bam']),
            mgd.TempOutputFile('normal.raw.bam'),
            mgd.TempOutputFile('tumour.unspiked.bam'),
            0.5,
        ),
    )

    workflow.commandline(
        name='simulate',
        args=(
            'destruct_bamextractsimreads',
            '-b', mgd.InputFile(args['bam']),
            '-r', mgd.InputFile(args['ref']),
            '-s', mgd.InputFile(os.path.join(args['outdir'], 'simulated.fasta')),
            '-f', mgd.TempInputObj('simulation.params').extract(lambda a: a['coverage_fraction']),
            '-1', mgd.OutputFile(os.path.join(args['outdir'], 'simulated.1.fastq')),
            '-2', mgd.OutputFile(os.path.join(args['outdir'], 'simulated.2.fastq')),
        ),
    )

    workflow.subworkflow(
        name='bwa_align',
        func=destruct.benchmark.align.bwa.workflow.bwa_align_workflow,
        args=(
            mgd.InputFile(args['ref']),
            mgd.InputFile(os.path.join(args['outdir'], 'simulated.1.fastq')),
            mgd.InputFile(os.path.join(args['outdir'], 'simulated.2.fastq')),
            mgd.TempOutputFile('simulated.unsorted.bam'),
        ),
    )

    workflow.transform(
        name='samtools_merge_sort_index',
        func=destruct.benchmark.destruct_test.samtools_merge_sort_index,
        args=(
            mgd.TempOutputFile('tumour_raw.bam'),
            mgd.TempInputFile('tumour.unspiked.bam'),
            mgd.TempInputFile('simulated.unsorted.bam'),
        ),
    )

    workflow.transform(
        name='samtools_reheader_tumour',
        func=destruct.benchmark.destruct_test.samtools_sample_reheader,
        args=(
            mgd.OutputFile(os.path.join(args['outdir'], 'tumour.bam')),
            mgd.TempInputFile('tumour_raw.bam'),
            'tumour',
        ),
    )

    workflow.transform(
        name='samtools_reheader_normal',
        func=destruct.benchmark.destruct_test.samtools_sample_reheader,
        args=(
            mgd.OutputFile(os.path.join(args['outdir'], 'normal.bam')),
            mgd.TempInputFile('normal.raw.bam'),
            'normal',
        ),
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('tool_defs', 'tool_name'),
        value=tool_defs,
    )

    workflow.subworkflow(
        name='run_tool',
        axes=('tool_name',),
        func=destruct.benchmark.destruct_test.create_tool_workflow,
        args=(
            mgd.TempInputObj('tool_defs', 'tool_name'),
            {
                'tumour': mgd.InputFile(os.path.join(args['outdir'], 'tumour.bam')),
                'normal': mgd.InputFile(os.path.join(args['outdir'], 'normal.bam')),
            },
            mgd.OutputFile(os.path.join(args['outdir'], 'results_{tool_name}.tsv'), 'tool_name'),
            mgd.TempSpace('tool_raw_data', 'tool_name'),
        ),
        kwargs={
            'control_id': 'normal',
        },
    )

    workflow.transform(
        name='plot',
        axes=('tool_name',),
        func=destruct.benchmark.destruct_test.create_roc_plot,
        args=(
            mgd.TempInputObj('simulation.params'),
            mgd.TempInputObj('tool_defs', 'tool_name'),
            mgd.InputFile(os.path.join(args['outdir'], 'simulated.tsv')),
            mgd.InputFile(os.path.join(args['outdir'], 'results_{tool_name}.tsv'), 'tool_name'),
            mgd.OutputFile(os.path.join(args['outdir'], 'annotated_{tool_name}.tsv'), 'tool_name'),
            mgd.OutputFile(os.path.join(args['outdir'], 'identified_{tool_name}.tsv'), 'tool_name'),
            mgd.OutputFile(os.path.join(args['outdir'], 'plots_{tool_name}.pdf'), 'tool_name'),
        ),
    )

    pyp.run(workflow)


