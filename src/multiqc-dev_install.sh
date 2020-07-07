#!/bin/bash

CONDA_ENV="/work/demassyie/bin/miniconda2/envs/SSDSnextflowPipeline"

eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV

cd ${NXF_PIPEDIR}/accessoryFiles/SSDS/MultiQC_SSDS_Rev1

python setup.py build
python setup.py install

conda deactivate
