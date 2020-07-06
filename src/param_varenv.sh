#!/bin/bash

CONDA_ENV="/work/demassyie/bin/miniconda2/envs/SSDSnextflowPipeline"
CONDA_NXF_GENOMES="/home/demassyie/work/SSDSnextflowPipeline/input_data/genomes"
CONDA_NXF_PIPEDIR="/home/demassyie/work/SSDSnextflowPipeline"
CONDA_PYTHONPATH="${CONDA_NXF_PIPEDIR}/accessoryFiles/SSDS/MultiQC_SSDS_Rev1/lib/python2.7"
CONDA_SCRATCH="/work/demassyie/scratch"

eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV
conda env config vars set NXF_GENOMES=${CONDA_NXF_GENOMES}
conda env config vars set NXF_PIPEDIR=${CONDA_NXF_PIPEDIR}
conda env config vars set PYTHONPATH=${CONDA_PYTHONPATH}
conda env config vars set SCRATCH=${CONDA_SCRATCH}
conda activate $CONDA_ENV

echo ${NXF_GENOMES}
echo ${NXF_PIPEDIR}
echo ${PYTHONPATH}
echo ${SCRATCH}

conda deactivate
echo ${NXF_GENOMES}
echo ${NXF_PIPEDIR}
echo ${PYTHONPATH}
echo ${SCRATCH}

