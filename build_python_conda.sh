
yum install gcc-c++ git make -y
conda config --set always_yes true
conda update conda
conda install python=3.9
conda config --add channels 'conda-forge'
conda config --add channels https://conda.anaconda.org/dranew
conda config --add channels https://conda.anaconda.org/shahcompbio
conda config --add channels 'bioconda'
conda install conda-build conda-verify anaconda-client
conda build conda/destruct --no-test
anaconda -t $CONDA_UPLOAD_TOKEN upload /usr/local/conda-bld/linux-64/destruct-*.tar.bz2


