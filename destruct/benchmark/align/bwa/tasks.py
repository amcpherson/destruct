import gzip
import itertools
import os

import pypeliner


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
            for name, seq, comment, qual in itertools.zip_longest(*[in_file]*4):
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


def merge_bam(input_bams, output_bam):
    if len(input_bams) >= 100:
        output_bam_tmp1 = output_bam + ".tmp1"
        output_bam_tmp2 = output_bam + ".tmp2"
        try: os.remove(output_bam_tmp1)
        except: pass
        try: os.remove(output_bam_tmp2)
        except: pass
        input_bams_1, input_bams_2 = split_dict(input_bams)
        merge_bam(input_bams_1, output_bam_tmp1)
        merge_bam(input_bams_2, output_bam_tmp2)
        merge_bam({'tmp1':output_bam_tmp1, 'tmp2':output_bam_tmp2}, output_bam)
        os.remove(output_bam_tmp1)
        os.remove(output_bam_tmp2)
    elif len(input_bams) == 1:
        pypeliner.commandline.execute('cp', list(input_bams.values())[0], output_bam)
    else:
        pypeliner.commandline.execute('samtools', 'merge', '-c', output_bam, *input_bams.values())
