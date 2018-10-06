FROM continuumio/miniconda
ARG app_version

RUN conda config --add channels https://conda.anaconda.org/dranew && conda config --add channels bioconda
RUN conda install destruct==$app_version
RUN mkdir /destruct_ref_data
RUN mkdir -p /root/.config/matplotlib
RUN echo "backend : Agg" > /root/.config/matplotlib/matplotlibrc
RUN destruct create_ref_data /destruct_ref_data

ENV NAME destruct

CMD ["destruct"]

