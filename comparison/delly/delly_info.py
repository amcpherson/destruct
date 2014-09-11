import os

import info
import utils


install_directory = os.path.join(info.install_directory, 'delly')

packages_directory = os.path.join(install_directory, 'packages')
bin_directory = os.path.join(install_directory, 'bin')

bamtools_root = os.path.join(packages_directory, 'bamtools')
seqtk_root = os.path.join(packages_directory, 'seqtk')
boost_root = install_directory

boost_url = 'http://downloads.sourceforge.net/project/boost/boost/1.56.0/boost_1_56_0.tar.gz'
boost_package = 'boost_1_56_0'
boost_tgz = boost_package+'.tar.gz'

delly_bin = os.path.join(bin_directory, 'delly')
