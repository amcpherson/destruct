#!/bin/bash

# Build and install executables
cd src
export CPLUS_INCLUDE_PATH=$PREFIX/include
scons install --prefix $PREFIX --boost_source $PREFIX/src/boost/
cd ..

# Build and install python package
$PYTHON setup.py install
