import csv
import sys
import logging
import os
import ConfigParser
import re
import itertools
import subprocess
import argparse
import string
import tarfile
import collections
import math
import pandas as pd
import numpy as np

import pygenes
import pypeliner
import pypeliner.workflow
import pypeliner.managed as mgd

import destruct.workflow


destruct_directory = os.environ.get('DESTRUCT_PACKAGE_DIRECTORY', None)
if destruct_directory is None:
    raise Exception('please set the $DESTRUCT_PACKAGE_DIRECTORY environment variable to the root of the destruct package')

data_directory = os.path.join(destruct_directory, 'data')
bin_directory = os.path.join(destruct_directory, 'bin')
default_config_filename = os.path.join(destruct_directory, 'defaultconfig.py')


__version__ = '0.1.0'

if __name__ == '__main__':

    argparser = argparse.ArgumentParser()

    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('--version', action='version', version=__version__)

    argparser.add_argument('ref_data_dir',
                           help='Reference dataset directory')

    argparser.add_argument('breakpoint_table',
                           help='Output table of breakpoint information in TSV format')

    argparser.add_argument('breakpoint_library_table',
                           help='Output table of library specific breakpoint information in TSV format')

    argparser.add_argument('--breakpoint_read_table', required=False, default=None,
                           help='Output table of breakpoint read information in TSV format')

    argparser.add_argument('--bam_files', nargs='+',
                           help='Input bam filenames')

    argparser.add_argument('--lib_ids', nargs='+',
                           help='Input ids for respective bam filenames')

    argparser.add_argument('--config', required=False,
                           help='Configuration Filename')

    args = vars(argparser.parse_args())

    if len(args['bam_files']) != len(args['lib_ids']):
        raise Exception('--lib_ids must correspond one to one with --bam_files')

    bam_filenames = dict(zip(args['lib_ids'], args['bam_files']))

    config = {'ref_data_directory':args['ref_data_dir'],
              'package_data_directory':data_directory}
    execfile(default_config_filename, {}, config)

    if args['config'] is not None:
        execfile(args['config'], {}, config)

    config.update(args)

    pyp = pypeliner.app.Pypeline([destruct], config)

    workflow = destruct.workflow.create_destruct_workflow(
        bam_filenames,
        args['breakpoint_table'],
        args['breakpoint_library_table'],
        config,
        bin_directory,
        breakpoint_read_table=args['breakpoint_read_table'],
    )

    pyp.run(workflow)


