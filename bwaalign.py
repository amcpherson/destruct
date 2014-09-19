import csv
import sys
import logging
import os
import ConfigParser
import itertools
import argparse
import string
import gzip

import pypeliner
import pypeliner.managed as mgd

if __name__ == '__main__':

    import bwaalign

    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pypeliner.app.add_arguments(argparser)
    argparser.add_argument('genome', help='Reference genome fasta')
    argparser.add_argument('fastq1', help='Fastq end 1 filename')
    argparser.add_argument('fastq2', help='Fastq end 2 filename')
    argparser.add_argument('output_bam', help='Output bam filename')
    argparser.add_argument('--reads_per_job', type=int, default=1000000, help='Reads per alignment job')
    argparser.add_argument('--config', help='Configuration filename')

    args = vars(argparser.parse_args())

    config = {}

    if args['config'] is not None:
        execfile(args['config'], {}, config)

    config.update(args)

    pyp = pypeliner.app.Pypeline([bwaalign], config)

    lowmem = {'mem':1}
    medmem = {'mem':8}
    himem = {'mem':32}

    if not os.path.exists(args['genome']+'.bwt'):
        raise Exception('No index for ' + args['genome'])

    pyp.sch.transform('split1', (), lowmem,
        bwaalign.split_fastq,
        None,
        mgd.InputFile(args['fastq1']),
        args['reads_per_job'],
        mgd.TempOutputFile('fastq1', 'byread'))

    pyp.sch.transform('split2', (), lowmem,
        bwaalign.split_fastq,
        None, 
        mgd.InputFile(args['fastq2']), 
        args['reads_per_job'], 
        mgd.TempOutputFile('fastq2', 'byread2'))

    pyp.sch.changeaxis('axis', (), 'fastq2', 'byread2', 'byread')

    pyp.sch.commandline('aln1', ('byread',), medmem,
        'bwa', 'aln',
        args['genome'],
        mgd.TempInputFile('fastq1', 'byread'),
        '>',
        mgd.TempOutputFile('sai1', 'byread'))

    pyp.sch.commandline('aln2', ('byread',), medmem,
        'bwa', 'aln',
        args['genome'],
        mgd.TempInputFile('fastq2', 'byread'),
        '>',
        mgd.TempOutputFile('sai2', 'byread'))

    pyp.sch.commandline('sampe', ('byread',), medmem,
        'bwa', 'sampe',
        args['genome'],
        mgd.TempInputFile('sai1', 'byread'),
        mgd.TempInputFile('sai2', 'byread'),
        mgd.TempInputFile('fastq1', 'byread'),
        mgd.TempInputFile('fastq2', 'byread'),
        '>',
        mgd.TempOutputFile('sam', 'byread'))

    pyp.sch.transform('cat', (), lowmem,
        bwaalign.cat_merge,
        None,
        mgd.TempInputFile('sam', 'byread'),
        mgd.TempOutputFile('sam'))

    pyp.sch.commandline('bam', (), lowmem,
        'grep', '-v', '^@',
        mgd.TempInputFile('sam'),
        '|',
        'samtools', 'view', '-bt',
        args['genome']+'.fai',
        '-',
        '>',
        mgd.OutputFile(args['output_bam']))

    pyp.run()

else:

    def split_fastq(in_filename, num_reads_per_file, out_filename_callback):
        if in_filename.endswith('.gz'):
            opener = gzip.open
        else:
            opener = open
        with opener(in_filename, 'r') as in_file:
            file_number = 0
            out_file = None
            out_file_read_count = None
            try:
                for name, seq, comment, qual in itertools.izip_longest(*[in_file]*4):
                    if out_file is None or out_file_read_count == num_reads_per_file:
                        if out_file is not None:
                            out_file.close()
                        out_file = open(out_filename_callback(file_number), 'w')
                        out_file_read_count = 0
                        file_number += 1
                    out_file.write(name)
                    out_file.write(seq)
                    out_file.write(comment)
                    out_file.write(qual)
                    out_file_read_count += 1
            finally:
                if out_file is not None:
                    out_file.close()

    def sort_bam(input_bam, output_bam):
        pypeliner.commandline.execute('samtools', 'sort', '-o', input_bam, output_bam+'.sortprefix', '>', output_bam)

    def split_dict(d):
        ditems = d.items()
        return dict(ditems[:len(ditems)/2]), dict(ditems[len(ditems)/2:])

    def cat_merge(in_filenames, out_filename):
        in_filename_args = list((a[1] for a in sorted(in_filenames.items())))
        cat_args = ['cat'] + in_filename_args + ['>', out_filename]
        pypeliner.commandline.execute(*cat_args)

    def create_bam(sam_filenames, bam_filename, samtools_bin, genome_fasta_fai):
        create_bam_args = ['cat'] + sam_filename_args + ['|', samtools_bin, 'view', '-bt', genome_fasta_fai, '-', '>', bam_filename]
        pypeliner.commandline.execute(*create_bam_args)

    def merge_bam(samtools_bin, input_bams, output_bam):
        if len(input_bams) >= 100:
            output_bam_tmp1 = output_bam + ".tmp1"
            output_bam_tmp2 = output_bam + ".tmp2"
            try: os.remove(output_bam_tmp1)
            except: pass
            try: os.remove(output_bam_tmp2)
            except: pass
            input_bams_1, input_bams_2 = split_dict(input_bams)
            merge_bam(samtools_bin, input_bams_1, output_bam_tmp1)
            merge_bam(samtools_bin, input_bams_2, output_bam_tmp2)
            merge_bam(samtools_bin, {'tmp1':output_bam_tmp1, 'tmp2':output_bam_tmp2}, output_bam)
            os.remove(output_bam_tmp1)
            os.remove(output_bam_tmp2)
        elif len(input_bams) == 1:
            pypeliner.commandline.execute('cp', input_bams.values()[0], output_bam)
        else:
            pypeliner.commandline.execute(samtools_bin, 'merge', output_bam, *input_bams.values())
