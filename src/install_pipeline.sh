#!/bin/bash

cd /home/${USER}/work/SSDSnextflowPipeline

eval "$(conda shell.bash hook)"

## Create conda environment to install nextflow
conda create -name nextflow-dev
conda activate nextflow-dev
conda install -c bioconda nextflow=20.04.1
conda deactivate 

## Create SSDSnextflowPipeline conda environment from yml file
eval "$(conda shell.bash hook)"
conda env create -f environment.yml
conda activate SSDSnextflowPipeline

## Install Multiqc_dev python library in the conda env
cd accessoryFiles/SSDS/MultiQC_SSDS_Rev1
python setup.py build
python setup.py install

conda deactivate
