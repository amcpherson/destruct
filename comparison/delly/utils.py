import os
import shutil
import subprocess


class Sentinal(object):
    def __init__(self, sentinal_filename):
        self.sentinal_filename = sentinal_filename
    def __call__(self, name):
        return Sentinal(self.sentinal_filename + name)
    @property
    def unfinished(self):
        if os.path.exists(self.sentinal_filename):
            print 'sentinal file ' + self.sentinal_filename + ' exists'
            return False
        return True
    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is None:
            with open(self.sentinal_filename, 'w') as sentinal_file:
                pass


def makedirs(directory):
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != 17:
            raise


def rmtree(directory):
    try:
        shutil.rmtree(directory)
    except OSError as e:
        if e.errno != 2:
            raise


def remove(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != 2:
            raise


def symlink(filename):
    link_name = os.path.join(os.getcwd(), os.path.basename(filename))
    remove(link_name)
    os.symlink(filename, link_name)


class CurrentDirectory(object):
    def __init__(self, directory):
        self.directory = directory
    def __enter__(self):
        self.prev_directory = os.getcwd()
        makedirs(self.directory)
        os.chdir(self.directory)
    def __exit__(self, *args):
        os.chdir(self.prev_directory)


def wget_file(url, filename):
    makedirs(os.path.dirname(filename))
    subprocess.check_call(['wget', '--no-check-certificate', url, '-O', filename])


def download_single_sra_dataset(sra_id, output_dir):
    url = 'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/{0}/{1}/{2}/{2}.sra'.format(sra_id[:3], sra_id[:6], sra_id)
    with CurrentDirectory(output_dir):
        subprocess.check_call(['wget', '-r', '-nH', '--cut-dirs=9', url])


def add_read_end(fastq_filename, read_end):
    temp_fastq_filename = fastq_filename + '.tmp'
    with open(fastq_filename, 'r') as f_in, open(temp_fastq_filename, 'w') as f_out:
        for line_number, line in enumerate(f_in):
            if line_number % 4 == 0:
                line = line.rstrip() + '/' + read_end + '\n'
            f_out.write(line)
    os.remove(fastq_filename)
    os.rename(temp_fastq_filename, fastq_filename)


def extract_single_sra_dataset(sra_id, output_dir):
    with CurrentDirectory(output_dir):
        subprocess.check_call(['fastq-dump', '-F', '--split-files', sra_id+'.sra'])
        for read_end in ('1', '2'):
            add_read_end(sra_id+'_'+read_end+'.fastq', read_end)

