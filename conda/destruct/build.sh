#!/bin/bash

# Build and install executables
cd src
scons install --prefix $PREFIX --boost_source $PREFIX/src/boost/
cd ..

# Build and install python package
python setup.py install
