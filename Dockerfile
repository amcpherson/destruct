FROM continuumio/miniconda3
ARG app_version

RUN apt-get update
RUN apt-get install g++ libncurses5 -y
RUN conda config --add channels https://conda.anaconda.org/dranew && conda config --add channels bioconda
RUN conda install destruct_utils==$app_version
RUN conda install openssl=1.0
RUN conda install bowtie dwgsim bwa samtools
RUN pip install destruct==$app_version
RUN mkdir -p /root/.config/matplotlib
RUN echo "backend : Agg" > /root/.config/matplotlib/matplotlibrc

ENV NAME destruct

ENTRYPOINT ["destruct"]

