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


destruct_directory = os.environ.get('DESTRUCT_PACKAGE_DIRECTORY', None)
if destruct_directory is None:
    raise Exception('please set the $DESTRUCT_PACKAGE_DIRECTORY environment variable to the root of the destruct package')

bin_directory = os.path.join(destruct_directory, 'bin')
default_config_filename = os.path.join(destruct_directory, 'defaultconfig.py')


if __name__ == '__main__':

    import destruct_test
    import destruct_bam_test
    import create_breakpoint_simulation

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)
    argparser.add_argument('simconfig',
                           help='Simulation configuration filename')

    argparser.add_argument('bam',
                           help='Source bam filename')

    argparser.add_argument('ref',
                           help='Reference genome for source bam')

    argparser.add_argument('installdir',
                           help='Tool installations directory')

    argparser.add_argument('outdir',
                           help='Output directory')

    argparser.add_argument('--config',
                           help='Configuration filename')

    args = vars(argparser.parse_args())

    config = {}

    if args['config'] is not None:
        execfile(args['config'], {}, config)

    config.update(args)

    pyp = pypeliner.app.Pypeline([destruct_test, destruct_bam_test, create_breakpoint_simulation], config)

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

    workflow.transform(
        name='create_sim',
        func=create_breakpoint_simulation.create_breakpoints,
        args=(
            mgd.TempInputObj('simulation.params'),
            mgd.InputFile(args['ref']),
            mgd.OutputFile(os.path.join(args['outdir'], 'simulated.fasta')),
            mgd.OutputFile(os.path.join(args['outdir'], 'simulated.tsv')),
        ),
    )

    workflow.transform(
        name='partition',
        func=destruct_test.partition_bam,
        args=(
            mgd.InputFile(args['bam']),
            mgd.OutputFile(os.path.join(args['outdir'], 'normal.bam')),
            mgd.TempOutputFile('tumour.unspiked.bam'),
            0.5,
        ),
    )        

    workflow.commandline(
        name='simulate',
        args=(
            os.path.join(bin_directory, 'bamextractsimreads'),
            '-b', mgd.InputFile(args['bam']),
            '-r', mgd.InputFile(args['ref']),
            '-s', mgd.InputFile(os.path.join(args['outdir'], 'simulated.fasta')),
            '-f', mgd.TempInputObj('simulation.params').extract(lambda a: a['coverage_fraction']),
            '-1', mgd.OutputFile(os.path.join(args['outdir'], 'simulated.1.fastq')),
            '-2', mgd.OutputFile(os.path.join(args['outdir'], 'simulated.2.fastq')),
        ),
    )        

    bwaalign_script = os.path.join(destruct_directory, 'scripts', 'bwaalign.py')

    workflow.commandline(
        name='bwa_align',
        args=(
            sys.executable,
            bwaalign_script,
            mgd.InputFile(args['ref']),
            mgd.InputFile(os.path.join(args['outdir'], 'simulated.1.fastq')),
            mgd.InputFile(os.path.join(args['outdir'], 'simulated.2.fastq')),
            mgd.TempOutputFile('simulated.unsorted.bam'),
            '--tmp', mgd.TempSpace('bwa_tmp'),
        ),
    )        

    workflow.transform(
        name='samtools_merge_sort_index',
        func=destruct_test.samtools_merge_sort_index,
        args=(
            mgd.OutputFile(os.path.join(args['outdir'], 'tumour.bam')),
            mgd.TempInputFile('tumour.unspiked.bam'),
            mgd.TempInputFile('simulated.unsorted.bam'),
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
            'control_id': 'normal',
            'tumour': mgd.InputFile(os.path.join(args['outdir'], 'tumour.bam')),
            'normal': mgd.InputFile(os.path.join(args['outdir'], 'normal.bam')),
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

