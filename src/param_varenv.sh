#!/bin/bash

CONDA_ENV="/work/demassyie/bin/miniconda2/envs/SSDSnextflowPipeline"
CONDA_NXF_GENOMES="/home/demassyie/work/SSDSnextflowPipeline/input_data/genomes"
CONDA_NXF_PIPEDIR="/home/demassyie/work/SSDSnextflowPipeline"
CONDA_ENV_PYTHON2="/work/demassyie/bin/miniconda2/envs/python2.7"
CONDA_SCRATCH="/work/demassyie/scratch"

eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV
conda env config vars set CONDA_ENV=${CONDA_ENV}
conda env config vars set NXF_GENOMES=${CONDA_NXF_GENOMES}
conda env config vars set NXF_PIPEDIR=${CONDA_NXF_PIPEDIR}
conda env config vars set PYTHON2=${CONDA_ENV_PYTHON2}
conda env config vars set SCRATCH=${CONDA_SCRATCH}
conda activate $CONDA_ENV
conda deactivate

conda activate $CONDA_ENV_PYTHON2
conda env config vars set NXF_PIPEDIR=${CONDA_NXF_PIPEDIR}
conda deactivate