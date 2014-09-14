import glob
import shutil
import os
import sys
import subprocess
import tarfile
import argparse

import utils


class DellyWrapper(object):

    def __init__(self, install_directory):

        self.install_directory = install_directory

        self.packages_directory = os.path.join(self.install_directory, 'packages')
        self.bin_directory = os.path.join(self.install_directory, 'bin')

        self.bamtools_root = os.path.join(self.packages_directory, 'bamtools')
        self.seqtk_root = os.path.join(self.packages_directory, 'seqtk')
        self.boost_root = self.install_directory

        self.boost_url = 'http://downloads.sourceforge.net/project/boost/boost/1.56.0/boost_1_56_0.tar.gz'
        self.boost_package = 'boost_1_56_0'
        self.boost_tgz = self.boost_package+'.tar.gz'

        self.delly_bin = os.path.join(self.bin_directory, 'delly')
        self.delly_excl_chrom = os.path.join(self.packages_directory, 'delly', 'human.hg19.excl.tsv')


    def install(self):

        Sentinal = utils.Sentinal(os.path.join(self.install_directory, 'sentinal_'))

        with Sentinal('install_bamtools') as sentinal:

            if sentinal.unfinished:

                with utils.CurrentDirectory(self.packages_directory):

                    utils.rmtree('bamtools')
                    subprocess.check_call('git clone https://github.com/pezmaster31/bamtools.git', shell=True)

                    with utils.CurrentDirectory('bamtools'):

                        utils.makedirs('build')

                        with utils.CurrentDirectory('build'):

                            subprocess.check_call('cmake ..', shell=True)
                            subprocess.check_call('make', shell=True)


        with Sentinal('install_kseq') as sentinal:

            if sentinal.unfinished:

                with utils.CurrentDirectory(self.packages_directory):

                    utils.rmtree('seqtk')
                    subprocess.check_call('git clone https://github.com/lh3/seqtk.git', shell=True)


        with Sentinal('install_boost') as sentinal:

            if sentinal.unfinished:

                with utils.CurrentDirectory(self.packages_directory):

                    subprocess.check_call('wget ' + self.boost_url, shell=True)
                    subprocess.check_call('tar -xzvf ' + self.boost_tgz, shell=True)

                    with utils.CurrentDirectory(self.boost_package):

                        subprocess.check_call('./bootstrap.sh', shell=True)

                        clean_env = 'INCLUDE_PATH= CPLUS_INCLUDE_PATH='

                        boost_build_cmd = clean_env + ' ./b2 install'
                        boost_build_cmd += ' link=static'
                        boost_build_cmd += ' --prefix=' + self.install_directory
                        boost_build_cmd += ' --with-iostreams'
                        boost_build_cmd += ' --with-filesystem'
                        boost_build_cmd += ' --with-program_options'
                        boost_build_cmd += ' --with-date_time'

                        subprocess.check_call(boost_build_cmd, shell=True)


        with Sentinal('download_delly') as sentinal:

            if sentinal.unfinished:

                with utils.CurrentDirectory(self.packages_directory):

                    utils.rmtree('delly')
                    subprocess.check_call('git clone https://github.com/tobiasrausch/delly.git', shell=True)

                    with utils.CurrentDirectory('delly'):

                        subprocess.check_call('git checkout v0.5.9', shell=True)

                        with open('Makefile.tmp', 'w') as f:
                            subprocess.check_call('sed s/-O9/-g\ -O3/g Makefile', shell=True, stdout=f)
                        os.rename('Makefile.tmp', 'Makefile')


        with Sentinal('install_delly') as sentinal:

            if sentinal.unfinished:

                with utils.CurrentDirectory(os.path.join(self.packages_directory, 'delly')):

                    make_cmd = 'make -B'
                    make_cmd += ' BAMTOOLS_ROOT=' + self.bamtools_root
                    make_cmd += ' SEQTK_ROOT=' + self.seqtk_root
                    make_cmd += ' BOOST_ROOT=' + self.boost_root
                    make_cmd += ' src/delly'

                    subprocess.check_call(make_cmd, shell=True)

            with utils.CurrentDirectory(self.bin_directory):

                utils.symlink(os.path.join(self.packages_directory, 'delly', 'src', 'delly'))


    def run_sv_type(self, sv_type, genome_fasta, bam_filenames, output_filename):

        with utils.SafeWriteFile(output_filename) as temp_output_filename:

            delly_cmd = list()
            delly_cmd += [self.delly_bin]
            delly_cmd += ['-t', sv_type]
            delly_cmd += ['-x', self.delly_excl_chrom]
            delly_cmd += ['-o', temp_output_filename]
            delly_cmd += ['-g', genome_fasta]
            delly_cmd += bam_filenames

            subprocess.check_call(delly_cmd)


    def convert_output(output_filename, ):
        pass


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('install_directory', help='Delly installation directory')
    args = parser.parse_args()

    delly = DellyWrapper(args.install_directory)

    delly.install()





