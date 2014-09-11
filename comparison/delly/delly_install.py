import glob
import shutil
import os
import sys
import subprocess
import tarfile

import utils
import info
import delly_info


Sentinal = utils.Sentinal(os.path.join(delly_info.install_directory, 'sentinal_'))


with Sentinal('install_bamtools') as sentinal:

    if sentinal.unfinished:

        with utils.CurrentDirectory(delly_info.packages_directory):

            utils.rmtree('bamtools')
            subprocess.check_call('git clone https://github.com/pezmaster31/bamtools.git', shell=True)

            with utils.CurrentDirectory('bamtools'):

                utils.makedirs('build')

                with utils.CurrentDirectory('build'):

                    subprocess.check_call('cmake ..', shell=True)
                    subprocess.check_call('make', shell=True)


with Sentinal('install_kseq') as sentinal:

    if sentinal.unfinished:

        with utils.CurrentDirectory(delly_info.packages_directory):

            utils.rmtree('seqtk')
            subprocess.check_call('git clone https://github.com/lh3/seqtk.git', shell=True)


with Sentinal('install_boost') as sentinal:

    if sentinal.unfinished:

        with utils.CurrentDirectory(delly_info.packages_directory):

            subprocess.check_call('wget ' + delly_info.boost_url, shell=True)
            subprocess.check_call('tar -xzvf ' + delly_info.boost_tgz, shell=True)

            with utils.CurrentDirectory(delly_info.boost_package):

                subprocess.check_call('./bootstrap.sh', shell=True)

                clean_env = 'INCLUDE_PATH= CPLUS_INCLUDE_PATH='

                boost_build_cmd = clean_env + ' ./b2 install'
                boost_build_cmd += ' link=static'
                boost_build_cmd += ' --prefix=' + delly_info.install_directory
                boost_build_cmd += ' --with-iostreams'
                boost_build_cmd += ' --with-filesystem'
                boost_build_cmd += ' --with-program_options'
                boost_build_cmd += ' --with-date_time'

                subprocess.check_call(boost_build_cmd, shell=True)


with Sentinal('install_delly') as sentinal:

    if sentinal.unfinished:

        with utils.CurrentDirectory(delly_info.packages_directory):

            utils.rmtree('delly')
            subprocess.check_call('git clone https://github.com/tobiasrausch/delly.git', shell=True)

            with utils.CurrentDirectory('delly'):

                subprocess.check_call('git checkout v0.5.9', shell=True)

                with open('Makefile.tmp', 'w') as f:
                    subprocess.check_call('sed s/-O9/-O3/g Makefile', shell=True, stdout=f)
                os.rename('Makefile.tmp', 'Makefile')

                make_cmd = 'make -B'
                make_cmd += ' BAMTOOLS_ROOT=' + delly_info.bamtools_root
                make_cmd += ' SEQTK_ROOT=' + delly_info.seqtk_root
                make_cmd += ' BOOST_ROOT=' + delly_info.boost_root
                make_cmd += ' src/delly'

                subprocess.check_call(make_cmd, shell=True)

        with utils.CurrentDirectory(delly_info.bin_directory):

            utils.symlink(os.path.join(delly_info.packages_directory, 'delly', 'src', 'delly'))



