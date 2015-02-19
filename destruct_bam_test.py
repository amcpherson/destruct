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
import pypeliner.managed as mgd

import wrappers
import utils.download

destruct_directory = os.path.abspath(os.path.dirname(__file__))

tools_directory = os.path.join(destruct_directory, 'tools')
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

    pyp.sch.transform('read_params', (), ctx,
        destruct_test.read_simulation_params,
        mgd.TempOutputObj('simulation.params'),
        mgd.InputFile(args['simconfig']))

    pyp.sch.transform('create_sim', (), ctx,
        create_breakpoint_simulation.create_breakpoints,
        None,
        mgd.TempInputObj('simulation.params'),
        mgd.InputFile(args['ref']),
        mgd.OutputFile(os.path.join(args['outdir'], 'simulated.fasta')),
        mgd.OutputFile(os.path.join(args['outdir'], 'simulated.tsv')))

    pyp.sch.transform('partition', (), ctx,
        destruct_test.partition_bam,
        None,
        mgd.InputFile(args['bam']),
        mgd.OutputFile(os.path.join(args['outdir'], 'normal.bam')),
        mgd.TempOutputFile('tumour.unspiked.bam'),
        0.5)

    pyp.sch.commandline('simulate', (), ctx,
        os.path.join(tools_directory, 'bamextractsimreads'),
        '-b', mgd.InputFile(args['bam']),
        '-r', mgd.InputFile(args['ref']),
        '-s', mgd.InputFile(os.path.join(args['outdir'], 'simulated.fasta')),
        '-f', mgd.TempInputObj('simulation.params').extract(lambda a: a['coverage_fraction']),
        '-1', mgd.OutputFile(os.path.join(args['outdir'], 'simulated.1.fastq')),
        '-2', mgd.OutputFile(os.path.join(args['outdir'], 'simulated.2.fastq')))

    bwaalign_script = os.path.join(destruct_directory, 'bwaalign.py')

    pyp.sch.commandline('bwa_align', (), ctx,
        sys.executable,
        bwaalign_script,
        mgd.InputFile(args['ref']),
        mgd.InputFile(os.path.join(args['outdir'], 'simulated.1.fastq')),
        mgd.InputFile(os.path.join(args['outdir'], 'simulated.2.fastq')),
        mgd.TempOutputFile('simulated.unsorted.bam'),
        '--tmp', mgd.TempFile('bwa_tmp'))

    pyp.sch.transform('samtools_merge_sort_index', (), ctx,
        destruct_test.samtools_merge_sort_index,
        None,
        mgd.OutputFile(os.path.join(args['outdir'], 'tumour.bam')),
        mgd.TempInputFile('tumour.unspiked.bam'),
        mgd.TempInputFile('simulated.unsorted.bam'))

    pyp.sch.transform('create_tool_wrappers', (), ctx,
        destruct_test.create_tool_wrappers,
        mgd.TempOutputObj('tool_wrapper', 'bytool'),
        args['installdir'])

    pyp.sch.transform('run_tool', ('bytool',), ctx,
        destruct_test.run_tool,
        None,
        mgd.TempInputObj('tool_wrapper', 'bytool'),
        mgd.TempFile('tool_tmp', 'bytool'),
        mgd.OutputFile(os.path.join(args['outdir'], 'results_{bytool}.tsv'), 'bytool'),
        control_id='normal',
        **{'tumour':mgd.InputFile(os.path.join(args['outdir'], 'tumour.bam')),
           'normal':mgd.InputFile(os.path.join(args['outdir'], 'normal.bam'))})


    pyp.sch.transform('plot', ('bytool',), ctx,
        destruct_test.create_roc_plot,
        None,
        mgd.TempInputObj('simulation.params'),
        mgd.TempInputObj('tool_wrapper', 'bytool'),
        mgd.InputFile(os.path.join(args['outdir'], 'simulated.tsv')),
        mgd.InputFile(os.path.join(args['outdir'], 'results_{bytool}.tsv'), 'bytool'),
        mgd.OutputFile(os.path.join(args['outdir'], 'annotated_{bytool}.tsv'), 'bytool'),
        mgd.OutputFile(os.path.join(args['outdir'], 'identified_{bytool}.tsv'), 'bytool'),
        mgd.OutputFile(os.path.join(args['outdir'], 'plots_{bytool}.pdf'), 'bytool'))

    pyp.run()


