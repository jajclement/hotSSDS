#!/bin/bash

cd /home/${USER}/work/ssdsnextflowpipeline

eval "$(conda shell.bash hook)"

## Create conda environment to install nextflow
conda create -name nextflow-dev
conda activate nextflow-dev
conda install -c bioconda nextflow=20.04.1
conda deactivate 

## Create ssdsnextflowpipeline conda environment from yml file
eval "$(conda shell.bash hook)"
conda env create -f environment.yml
conda activate ssdsnextflowpipeline

## Install Multiqc_dev python library in the conda env
cd accessoryFiles/SSDS/MultiQC_SSDS_Rev1
python setup.py build
python setup.py install

conda deactivate
