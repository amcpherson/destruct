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

    cfg = pypeliner.easypypeliner.Config(vars(argparser.parse_args()))
    pyp = pypeliner.easypypeliner.EasyPypeliner([], cfg)

    class GenesisQsubJob(pypeliner.execqueue.QsubJob):
        def create_submit_command(self, ctx, name, script_filename, qsub_bin, native_spec, stdout_filename, stderr_filename):
            if 'txfer' in ctx:
                qsub = 'qsub -sync y -b y -P transfer -q thosts.q -S /bin/bash'
                qsub += ' -N ' + name
                qsub += ' ' + script_filename
                return ['ssh', 'apollo', qsub]
            else:
                return super(GenesisQsubJob, self).create_submit_command(ctx, name, script_filename, qsub_bin, native_spec, stdout_filename, stderr_filename)

    class GenesisQsubJobQueue(pypeliner.execqueue.QsubJobQueue):
        def create(self, ctx, job):
            return GenesisQsubJob(ctx, job, self.temps_dir, self.modules, self.qsub_bin, self.native_spec)

    pyp.exec_queue = GenesisQsubJobQueue(pyp.exc_dir, [destruct, destruct_genesis], cfg.nativespec)

    txfer = {'txfer':True}

    pyp.sch.transform('readlibs', (), destruct.lowmem, destruct_genesis.read_libraries, pyp.sch.oobj('bamremote', ('bylibrary',)), pyp.sch.input(cfg.libraries))
    pyp.sch.commandline('copybam', ('bylibrary',), txfer, 'cp', pyp.sch.iobj('bamremote', ('bylibrary',)), pyp.sch.ofile('bam', ('bylibrary',)))
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

