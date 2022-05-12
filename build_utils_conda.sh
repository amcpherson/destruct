
yum install gcc-c++ libquadmath-devel git -y
conda config --set always_yes true
conda config --add channels https://conda.anaconda.org/dranew
conda config --add channels https://conda.anaconda.org/shahcompbio
conda config --add channels 'bioconda'
conda install conda-build
conda build conda/destruct_utils
conda install anaconda-client==1.6.7
anaconda -t $CONDA_UPLOAD_TOKEN upload /usr/local/conda-bld/linux-64/destruct*.tar.bz2
