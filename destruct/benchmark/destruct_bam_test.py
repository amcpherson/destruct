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

    argparser.add_argument('sim_config',
                           help='Simulation configuration filename')

    argparser.add_argument('ref_data_dir',
                           help='Reference genomes directory')

    argparser.add_argument('bam',
                           help='Source bam filename')

    argparser.add_argument('ref',
                           help='Reference genome for source bam')

    argparser.add_argument('results_dir',
                           help='Output directory')

    argparser.add_argument('tool_defs',
                           help='Tool Definition Filename')

    args = vars(argparser.parse_args())

    try:
        os.makedirs(args['results_dir'])
    except OSError:
        pass

    args['tmpdir'] = os.path.join(args['results_dir'], 'tmp')

    pyp = pypeliner.app.Pypeline(config=args)

    yaml_text = open(args['tool_defs']).read().format(ref_data_dir=args['ref_data_dir'])
    tool_defs = yaml.load(yaml_text)

    sim_config = yaml.load(open(args['sim_config']))

    ctx = {'mem':4}

    workflow = pypeliner.workflow.Workflow(default_ctx=ctx)

    workflow.setobj(mgd.TempOutputObj('simulation.params'), sim_config['simulation'])

    if sim_config['reference'].get('chromosomes', None) is not None:
        genome_fasta = os.path.join(args['results_dir'], 'genome.fa')
        source_bam = os.path.join(args['results_dir'], 'source.bam')

        workflow.setobj(mgd.TempOutputObj('chromosomes'), sim_config['reference']['chromosomes'])

        workflow.transform(
            name='create_ref_bam',
            func=destruct.benchmark.destruct_test.create_ref_bam,
            args=(
                mgd.InputFile(args['ref']),
                mgd.InputFile(args['bam']),
                mgd.OutputFile(genome_fasta),
                mgd.OutputFile(source_bam),
                mgd.TempInputObj('chromosomes'),
            ),
        )
    
    else:
        genome_fasta = args['ref']
        source_bam = args['bam']

    workflow.transform(
        name='create_sim',
        func=destruct.benchmark.create_breakpoint_simulation.create_breakpoints,
        args=(
            mgd.TempInputObj('simulation.params'),
            mgd.InputFile(genome_fasta),
            mgd.OutputFile(os.path.join(args['results_dir'], 'simulated.fa')),
            mgd.OutputFile(os.path.join(args['results_dir'], 'simulated.tsv')),
        ),
    )

    workflow.transform(
        name='partition',
        func=destruct.benchmark.destruct_test.partition_bam,
        args=(
            mgd.InputFile(source_bam),
            mgd.TempOutputFile('normal.raw.bam'),
            mgd.TempOutputFile('tumour.unspiked.bam'),
            0.5,
        ),
    )

    workflow.commandline(
        name='simulate',
        args=(
            'destruct_bamextractsimreads',
            '-b', mgd.InputFile(source_bam),
            '-r', mgd.InputFile(genome_fasta),
            '-s', mgd.InputFile(os.path.join(args['results_dir'], 'simulated.fa')),
            '-f', mgd.TempInputObj('simulation.params').extract(lambda a: a['coverage_fraction']),
            '-1', mgd.OutputFile(os.path.join(args['results_dir'], 'simulated.1.fastq')),
            '-2', mgd.OutputFile(os.path.join(args['results_dir'], 'simulated.2.fastq')),
        ),
    )

    workflow.subworkflow(
        name='bwa_align',
        func=destruct.benchmark.align.bwa.workflow.bwa_align_workflow,
        args=(
            mgd.InputFile(genome_fasta),
            mgd.InputFile(os.path.join(args['results_dir'], 'simulated.1.fastq')),
            mgd.InputFile(os.path.join(args['results_dir'], 'simulated.2.fastq')),
            mgd.TempOutputFile('simulated.unsorted.bam'),
        ),
        kwargs={
            'read_group_str': '@RG\\tID:B',
        },
    )

    workflow.transform(
        name='samtools_merge_sort_index',
        func=destruct.benchmark.destruct_test.samtools_merge_sort_index,
        args=(
            mgd.TempOutputFile('tumour.raw.bam'),
            mgd.TempInputFile('tumour.unspiked.bam'),
            mgd.TempInputFile('simulated.unsorted.bam'),
        ),
    )

    workflow.transform(
        name='samtools_reheader_tumour',
        func=destruct.benchmark.destruct_test.samtools_sample_reheader,
        args=(
            mgd.OutputFile(os.path.join(args['results_dir'], 'tumour.bam')),
            mgd.TempInputFile('tumour.raw.bam'),
            'tumour',
        ),
    )

    workflow.transform(
        name='samtools_reheader_normal',
        func=destruct.benchmark.destruct_test.samtools_sample_reheader,
        args=(
            mgd.OutputFile(os.path.join(args['results_dir'], 'normal.bam')),
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
                'tumour': mgd.InputFile(os.path.join(args['results_dir'], 'tumour.bam')),
                'normal': mgd.InputFile(os.path.join(args['results_dir'], 'normal.bam')),
            },
            mgd.OutputFile(os.path.join(args['results_dir'], 'results_{tool_name}.tsv'), 'tool_name'),
            mgd.TempSpace('tool_raw_data', 'tool_name', cleanup='none'),
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
            mgd.InputFile(genome_fasta),
            mgd.InputFile(os.path.join(args['results_dir'], 'simulated.tsv')),
            mgd.InputFile(os.path.join(args['results_dir'], 'results_{tool_name}.tsv'), 'tool_name'),
            mgd.OutputFile(os.path.join(args['results_dir'], 'annotated_{tool_name}.tsv'), 'tool_name'),
            mgd.OutputFile(os.path.join(args['results_dir'], 'identified_{tool_name}.tsv'), 'tool_name'),
            mgd.OutputFile(os.path.join(args['results_dir'], 'plots_{tool_name}.pdf'), 'tool_name'),
        ),
    )

    pyp.run(workflow)


