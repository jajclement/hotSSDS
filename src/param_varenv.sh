#!/bin/bash

CONDA_ENV="/work/demassyie/bin/miniconda2/envs/SSDSnextflowPipeline"
CONDA_NXF_GENOMES="/home/demassyie/work/SSDSnextflowPipeline/input_data/genomes"
CONDA_NXF_PIPEDIR="/home/demassyie/work/SSDSnextflowPipeline"
CONDA_ENV_PYTHON2="/work/demassyie/bin/miniconda2/envs/python2.7"
CONDA_SCRATCH="/work/demassyie/scratch"
CONDA_TOWER_ACCESS_TOKEN="db33ecd8b1cc2e1b39972b4d82614d379b84f878"
CONDA_NXF_VER="20.06.0-edge"

eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV
conda env config vars set CONDA_ENV=${CONDA_ENV}
conda env config vars set NXF_GENOMES=${CONDA_NXF_GENOMES}
conda env config vars set NXF_PIPEDIR=${CONDA_NXF_PIPEDIR}
conda env config vars set PYTHON2=${CONDA_ENV_PYTHON2}
conda env config vars set SCRATCH=${CONDA_SCRATCH}
conda env config vars set TOWER_ACCESS_TOKEN=${CONDA_TOWER_ACCESS_TOKEN}
conda env config vars set NXF_VER=${CONDA_NXF_VER}
conda activate $CONDA_ENV
conda deactivate

conda activate $CONDA_ENV_PYTHON2
conda env config vars set NXF_PIPEDIR=${CONDA_NXF_PIPEDIR}
conda deactivate
