import fnmatch
import os
import sys
import platform
import argparse

from distutils.core import setup
from distutils.extension import Extension

argparser = argparse.ArgumentParser(add_help=False)
argparser.add_argument('--boost_dir', help='boost source directory', required=True)
args, unknown = argparser.parse_known_args()
sys.argv = [sys.argv[0]] + unknown

def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename

boost_python_sources = list(find_files(os.path.join(args.boost_dir, 'libs/python/src'), '*.cpp'))

assert(len(boost_python_sources) > 0)

extra_compile_args = ['-Wno-unused-variable']
if platform.platform().startswith('Darwin'):
    extra_compile_args.append('-Wno-unneeded-internal-declaration')
    extra_compile_args.append('-Wno-unused-private-field')

extra_link_args = []
if platform.platform().startswith('Darwin'):
    extra_link_args.append('-Wl,-no_compact_unwind')

ext_modules = [Extension('pybreakcopies', 
                         ['pybreakcopies.cpp'] + boost_python_sources,
                         language='c++',
                         define_macros=[('PYTHON_PACKAGE', None)],
                         libraries=['ipopt', 'blas', 'lapack', 'gfortran', 'coinmumps', 'coinmetis'],
                         extra_compile_args=extra_compile_args,
                         extra_link_args=extra_link_args
                        )]

setup(name='pybreakcopies', ext_modules=ext_modules)

