import csv
import sys
import logging
import os
import ConfigParser
import re
import itertools
import collections
import subprocess
import argparse
import string

import pypeliner
import pypeliner.workflow
import pypeliner.managed as mgd

import wrappers
import utils.download

destruct_directory = os.path.abspath(os.path.dirname(__file__))

default_config_filename = os.path.join(destruct_directory, 'defaultconfig.py')


if __name__ == '__main__':

    import destruct_test
    import destruct_sim_test
    import create_breakpoint_simulation

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('simconfig',
                           help='Simulation configuration filename')

    argparser.add_argument('installdir',
                           help='Tool installations directory')

    argparser.add_argument('outdir',
                           help='Output directory')

    argparser.add_argument('--config',
                           help='Configuration filename')

    argparser.add_argument('--chromosomes', nargs='*', type=str, default=['20'],
                           help='Reference chromosomes')

    argparser.add_argument('--include_nonchromosomal',  action='store_true',
                           help='Include non chromosomal reference sequences')

    args = vars(argparser.parse_args())

    config = {}

    if args['config'] is not None:
        execfile(args['config'], {}, config)

    config.update(args)

    pyp = pypeliner.app.Pypeline([destruct_test, destruct_sim_test, create_breakpoint_simulation], config)

    try:
        os.makedirs(args['outdir'])
    except OSError:
        pass

    ctx = {'mem':4}

    workflow = pypeliner.workflow.Workflow(default_ctx=ctx)

    workflow.transform(
        name='read_params',
        func=destruct_test.read_simulation_params,
        ret=mgd.TempOutputObj('simulation.params'),
        args=(mgd.InputFile(args['simconfig']),),
    )

    workflow.setobj(mgd.TempOutputObj('chromosomes'), args['chromosomes'])
    workflow.setobj(mgd.TempOutputObj('include_nonchromosomal'), args['include_nonchromosomal'])

    workflow.transform(
        name='create_genome',
        func=destruct_test.create_genome,
        args=(
            mgd.TempInputObj('chromosomes'),
            mgd.TempInputObj('include_nonchromosomal'),
            mgd.OutputFile(os.path.join(args['outdir'], 'genome.fasta')),
        ),
    )

    workflow.transform(
        name='create_sim',
        func=create_breakpoint_simulation.create,
        args=(
            mgd.TempInputObj('simulation.params'),
            mgd.InputFile(os.path.join(args['outdir'], 'genome.fasta')),
            mgd.OutputFile(os.path.join(args['outdir'], 'simulated.fasta')),
            mgd.OutputFile(os.path.join(args['outdir'], 'simulated.tsv')),
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
            '>', mgd.OutputFile(os.path.join(args['outdir'], 'simulated.1.fastq')),
        ),
    )
    
    workflow.commandline(
        name='cat2',
        args=(
            'cat',
            mgd.TempInputFile('concordant.2.fastq'),
            mgd.TempInputFile('discordant.2.fastq'),
            '>', mgd.OutputFile(os.path.join(args['outdir'], 'simulated.2.fastq')),
        ),
    )

    bwaalign_script = os.path.join(destruct_directory, 'bwaalign.py')

    workflow.commandline(
        name='bwa_align',
        args=(
            sys.executable,
            bwaalign_script,
            mgd.InputFile(os.path.join(args['outdir'], 'genome.fasta')),
            mgd.InputFile(os.path.join(args['outdir'], 'simulated.1.fastq')),
            mgd.InputFile(os.path.join(args['outdir'], 'simulated.2.fastq')),
            mgd.TempOutputFile('simulated.unsorted.bam'),
            '--tmp', mgd.TempSpace('bwa_tmp'),
        ),
    )

    workflow.transform(
        name='samtools_sort_index',
        func=destruct_test.samtools_sort_index,
        args=(
            mgd.TempInputFile('simulated.unsorted.bam'),
            mgd.OutputFile(os.path.join(args['outdir'], 'simulated.bam')),
        ),
    )

    workflow.transform(
        name='create_tool_wrappers',
        func=destruct_test.create_tool_wrappers,
        ret=mgd.TempOutputObj('tool_wrapper', 'bytool'),
        args=(args['installdir'],),
    )

    workflow.transform(
        name='run_tool',
        axes=('bytool',),
        func=destruct_test.run_tool,
        args=(
            mgd.TempInputObj('tool_wrapper', 'bytool'),
            mgd.TempSpace('tool_tmp', 'bytool'),
            mgd.OutputFile(os.path.join(args['outdir'], 'results_{bytool}.tsv'), 'bytool'),
        ),
        kwargs={
            'simulated': mgd.InputFile(os.path.join(args['outdir'], 'simulated.bam')),
        },
    )

    workflow.transform(
        name='plot',
        axes=('bytool',),
        func=destruct_test.create_roc_plot,
        args=(
            mgd.TempInputObj('simulation.params'),
            mgd.TempInputObj('tool_wrapper', 'bytool'),
            mgd.InputFile(os.path.join(args['outdir'], 'simulated.tsv')),
            mgd.InputFile(os.path.join(args['outdir'], 'results_{bytool}.tsv'), 'bytool'),
            mgd.OutputFile(os.path.join(args['outdir'], 'annotated_{bytool}.tsv'), 'bytool'),
            mgd.OutputFile(os.path.join(args['outdir'], 'identified_{bytool}.tsv'), 'bytool'),
            mgd.OutputFile(os.path.join(args['outdir'], 'plots_{bytool}.pdf'), 'bytool'),
        ),
    )

    pyp.run(workflow)

