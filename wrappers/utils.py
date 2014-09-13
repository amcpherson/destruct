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


class SafeWriteFile(object):
    def __init__(self, filename):
        self.filename = filename
        self.temp_filename = filename + '.tmp'
    def __enter__(self):
        return self.temp_filename
    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is None:
            os.rename(self.temp_filename, self.filename)

