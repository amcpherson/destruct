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

    argparser.add_argument('results_dir',
                           help='Output directory')

    argparser.add_argument('tool_defs',
                           help='Tool Definition Filename')

    argparser.add_argument('--chromosomes', nargs='*', type=str, default=['20'],
                           help='Reference chromosomes')

    argparser.add_argument('--include_nonchromosomal',  action='store_true',
                           help='Include non chromosomal reference sequences')

    args = vars(argparser.parse_args())

    try:
        os.makedirs(args['results_dir'])
    except OSError:
        pass

    args['tmpdir'] = os.path.join(args['results_dir'], 'tmp')

    pyp = pypeliner.app.Pypeline(config=args)

    yaml_text = open(args['tool_defs']).read().format(ref_data_dir=args['ref_data_dir'])
    tool_defs = yaml.load(yaml_text)

    ctx = {'mem':4}

    workflow = pypeliner.workflow.Workflow(default_ctx=ctx)

    workflow.transform(
        name='read_params',
        func=destruct.benchmark.destruct_test.read_simulation_params,
        ret=mgd.TempOutputObj('simulation.params'),
        args=(mgd.InputFile(args['sim_config']),),
    )

    workflow.setobj(mgd.TempOutputObj('chromosomes'), args['chromosomes'])
    workflow.setobj(mgd.TempOutputObj('include_nonchromosomal'), args['include_nonchromosomal'])

    workflow.transform(
        name='create_genome',
        func=destruct.benchmark.destruct_test.create_genome,
        args=(
            mgd.TempInputObj('chromosomes'),
            mgd.TempInputObj('include_nonchromosomal'),
            mgd.OutputFile(os.path.join(args['results_dir'], 'genome.fasta')),
        ),
    )

    workflow.transform(
        name='create_sim',
        func=destruct.benchmark.create_breakpoint_simulation.create,
        args=(
            mgd.TempInputObj('simulation.params'),
            mgd.InputFile(os.path.join(args['results_dir'], 'genome.fasta')),
            mgd.OutputFile(os.path.join(args['results_dir'], 'simulated.fasta')),
            mgd.OutputFile(os.path.join(args['results_dir'], 'simulated.tsv')),
            mgd.TempOutputFile('concordant.1.fastq'),
            mgd.TempOutputFile('concordant.2.fastq'),
            mgd.TempOutputFile('discordant.1.fastq'),
            mgd.TempOutputFile('discordant.2.fastq'),
        ),
    )

    workflow.commandline(
        name='cat1',
        args=(
            'cat',
            mgd.TempInputFile('concordant.1.fastq'),
            mgd.TempInputFile('discordant.1.fastq'),
            '>', mgd.OutputFile(os.path.join(args['results_dir'], 'simulated.1.fastq')),
        ),
    )
    
    workflow.commandline(
        name='cat2',
        args=(
            'cat',
            mgd.TempInputFile('concordant.2.fastq'),
            mgd.TempInputFile('discordant.2.fastq'),
            '>', mgd.OutputFile(os.path.join(args['results_dir'], 'simulated.2.fastq')),
        ),
    )

    workflow.subworkflow(
        name='bwa_align',
        func=destruct.benchmark.align.bwa.workflow.bwa_align_workflow,
        args=(
            mgd.InputFile(os.path.join(args['results_dir'], 'genome.fasta')),
            mgd.InputFile(os.path.join(args['results_dir'], 'simulated.1.fastq')),
            mgd.InputFile(os.path.join(args['results_dir'], 'simulated.2.fastq')),
            mgd.TempOutputFile('simulated.unsorted.bam'),
        ),
    )

    workflow.transform(
        name='samtools_sort_index',
        func=destruct.benchmark.destruct_test.samtools_sort_index,
        args=(
            mgd.TempInputFile('simulated.unsorted.bam'),
            mgd.OutputFile(os.path.join(args['results_dir'], 'simulated.bam')),
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
            {'simulated': mgd.InputFile(os.path.join(args['results_dir'], 'simulated.bam')),},
            mgd.OutputFile(os.path.join(args['results_dir'], 'results_{tool_name}.tsv'), 'tool_name'),
            mgd.TempSpace('tool_raw_data', 'tool_name'),
        ),
    )

    workflow.transform(
        name='plot',
        axes=('tool_name',),
        func=destruct.benchmark.destruct_test.create_roc_plot,
        args=(
            mgd.TempInputObj('simulation.params'),
            mgd.TempInputObj('tool_defs', 'tool_name'),
            mgd.InputFile(os.path.join(args['results_dir'], 'simulated.tsv')),
            mgd.InputFile(os.path.join(args['results_dir'], 'results_{tool_name}.tsv'), 'tool_name'),
            mgd.OutputFile(os.path.join(args['results_dir'], 'annotated_{tool_name}.tsv'), 'tool_name'),
            mgd.OutputFile(os.path.join(args['results_dir'], 'identified_{tool_name}.tsv'), 'tool_name'),
            mgd.OutputFile(os.path.join(args['results_dir'], 'plots_{tool_name}.pdf'), 'tool_name'),
        ),
    )

    pyp.run(workflow)


