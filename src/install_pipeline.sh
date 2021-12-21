#!/bin/bash

#Default parameters
PIPELINE_DIRECTORY="/home/${USER}/work/ssdsnextflowpipeline"
MAIN_ENV_NAME="nextflow21"
NXF_V="21.10.0"

#Read arguments from command line
while getopts hf:p:n:v: flag
do
    case "${flag}" in
	h) echo "Installation script for nextflow conda environment" ; echo "Usage: bash `basename $0` [options]"; echo "Options : "; echo "-h display help message"; echo "-v Version of Nextflow to install in conda environment (default : ${NXF_V})"; echo "-p Path to ssdsnextflowpipeline directory (default : /home/${USER}/work/ssdsnextflowpipeline)"; echo "-n Name for Nextflow conda environment (default : ${MAIN_ENV_NAME})"; exit 0;;
        p) PIPELINE_DIRECTORY=${OPTARG};if [ ! -d ${PIPELINE_DIRECTORY} ]; then echo "Directory ${PIPELINE_DIRECTORY} not found!" ; exit 0; fi;;
        n) MAIN_ENV_NAME=${OPTARG};;
	v) NXF_V=${OPTARG};;
    esac
done

#Go to pipeline directory
cd ${PIPELINE_DIRECTORY}

eval "$(conda shell.bash hook)"

## Create conda environment to install nextflow
conda create -y -n ${MAIN_ENV_NAME}
conda activate ${MAIN_ENV_NAME}
conda install -c bioconda nextflow=${NXF_V} --yes
conda deactivate 

echo "Installation of ${MAIN_ENV_NAME} conda environment with Nextflow version ${NXF_V} complete."


