#!/bin/bash

# Build and install executables
cd src
export CPLUS_INCLUDE_PATH=$PREFIX/include
export LIBRARY_PATH=$PREFIX/lib
scons install --prefix $PREFIX
cd ..

