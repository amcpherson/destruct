#!/bin/bash

# Build executables
cd src
scons install
cd ..

# Environment variable for location of executables
export DESTRUCT_PACKAGE_DIRECTORY=`pwd`/

python setup.py install
