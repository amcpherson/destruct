docker run -v ${PWD}:/repo -w /repo -e CONDA_UPLOAD_TOKEN=$CONDA_UPLOAD_TOKEN -it conda/miniconda2-centos6 bash -e build_utils_conda.sh
