docker run -v ${PWD}:/repo -w /repo -e CONDA_UPLOAD_TOKEN=$CONDA_UPLOAD_TOKEN -it conda/miniconda3-centos7 bash -e build_python_conda.sh
