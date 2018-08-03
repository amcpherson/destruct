import os
import argparse
import yaml

import pypeliner
import pypeliner.workflow
import pypeliner.managed as mgd

import destruct.benchmark.align.bwa.workflow
import destruct.benchmark.destruct_test
import destruct.benchmark.create_breakpoint_simulation
import destruct.benchmark.generate_bam


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
    workflow.setobj(mgd.TempOutputObj('chromosomes'), sim_config['reference']['chromosomes'])
    workflow.setobj(mgd.TempOutputObj('include_nonchromosomal'), sim_config['reference']['include_nonchromosomal'])

    workflow.subworkflow(
        name='generate_bam_workflow',
        func=destruct.benchmark.generate_bam.generate_bam,
        args=(
            mgd.TempInputObj('simulation.params'),
            mgd.TempInputObj('chromosomes'),
            mgd.TempInputObj('include_nonchromosomal'),
            mgd.OutputFile(os.path.join(args['results_dir'], 'simulated.bam')),
            mgd.OutputFile(os.path.join(args['results_dir'], 'genome.fasta')),
            mgd.OutputFile(os.path.join(args['results_dir'], 'simulated.tsv')),
            os.path.join(args['results_dir'], 'raw'),
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
            mgd.TempSpace('tool_raw_data', 'tool_name', cleanup='none'),
        ),
    )

    workflow.transform(
        name='plot',
        axes=('tool_name',),
        func=destruct.benchmark.destruct_test.create_roc_plot,
        args=(
            mgd.TempInputObj('simulation.params'),
            mgd.TempInputObj('tool_defs', 'tool_name'),
            mgd.InputFile(os.path.join(args['results_dir'], 'genome.fasta')),
            mgd.InputFile(os.path.join(args['results_dir'], 'simulated.tsv')),
            mgd.InputFile(os.path.join(args['results_dir'], 'results_{tool_name}.tsv'), 'tool_name'),
            mgd.OutputFile(os.path.join(args['results_dir'], 'annotated_{tool_name}.tsv'), 'tool_name'),
            mgd.OutputFile(os.path.join(args['results_dir'], 'identified_{tool_name}.tsv'), 'tool_name'),
            mgd.OutputFile(os.path.join(args['results_dir'], 'plots_{tool_name}.pdf'), 'tool_name'),
        ),
    )

    pyp.run(workflow)
