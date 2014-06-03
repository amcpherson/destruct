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
from collections import *

import pygenes
import pypeliner

if __name__ == '__main__':

    import destruct
    import destruct_genesis

    argparser = argparse.ArgumentParser()
    pypeliner.easypypeliner.add_arguments(argparser)
    argparser.add_argument('libraries', help='Libraries list filename')
    argparser.add_argument('breakpoints', help='Breakpoints table filename')
    argparser.add_argument('breakreads', help='Breakpoint reads table filename')
    argparser.add_argument('plots_tar', help='Diagnostic plots tar filename')

    cfg = pypeliner.easypypeliner.Config(vars(argparser.parse_args()))
    pyp = pypeliner.easypypeliner.EasyPypeliner([destruct, destruct_genesis], cfg)

    txfer = {'local':True}

    pyp.sch.transform('readlibs', (), destruct.lowmem, destruct_genesis.read_libraries, pyp.sch.oobj('bamremote', ('bylibrary',)), pyp.sch.input(cfg.libraries))
    pyp.sch.transform('copybam', ('bylibrary',), txfer, destruct_genesis.copy_bam, None, pyp.sch.iobj('bamremote', ('bylibrary',)), pyp.sch.ofile('bam', ('bylibrary',)))
    destruct.multilib_predict_breakpoints(pyp.sch, cfg, pyp.sch.ifile('bam', ('bylibrary',)), pyp.sch.output(cfg.breakpoints), pyp.sch.output(cfg.breakreads))

    pyp.run()

else:

    def read_libraries(libraries_filename):
        libraries = dict()
        with open(libraries_filename, 'r') as libraries_file:
            for row in csv.reader(libraries_file, delimiter='\t'):
                lib_name = row[0]
                lib_bam = row[1]
                libraries[lib_name] = lib_bam
        return libraries

    def copy_bam(remote_path, local_path):
        qsub = 'qsub -sync y -b y -P transfer -q thosts.q -S /bin/bash'
        qsub += ' -N copybam'
        qsub += ' -b yes'
        qsub += ' cp ' + remote_path + ' ' + local_path
        pypeliner.commandline.execute(*['ssh', 'apollo', qsub])

