
yum install gcc-c++ git -y
conda config --set always_yes true
conda config --add channels https://conda.anaconda.org/dranew
conda config --add channels https://conda.anaconda.org/shahcompbio
conda config --add channels 'bioconda'
conda install -y scons boost_lib==1.60.0 tclap==1.2.1 bzip2

cd src
export CPLUS_INCLUDE_PATH=$PREFIX/include
export LIBRARY_PATH=$PREFIX/lib
scons build
cd ..

export LD_LIBRARY_PATH=/usr/local/lib/

