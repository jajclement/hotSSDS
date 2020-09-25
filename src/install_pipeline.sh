#!/bin/bash

#Default parameters
PIPELINE_DIRECTORY="/home/${USER}/work/ssdsnextflowpipeline"
MULTIQC_DEV_DIRECTORY="${PIPELINE_DIRECTORY}/accessoryFiles/SSDS/MultiQC_SSDS_Rev1"
MAIN_ENV_NAME="nextflow_dev"
MULTIQC_DEV_NAME="multiqc_dev"
FLAG=1

#Read arguments from command line
while getopts hf:p:m:e:c: flag
do
    case "${flag}" in
	h) echo "Installation script for nextflow and multiqc-dev conda environments" ; echo "Usage: bash `basename $0` [options]"; echo "Options : "; echo "-h display help message"; echo "-f value should be 0 (install none of the conda environments) ; 1 (install only Nextflow conda environment) ; 2 (install both nextflow and multiqc-dev environments) or 3 (install only multiqc-dev conda environment) (default 1)"; echo "-p Path to ssdsnextflowpipeline directory (default : /home/${USER}/work/ssdsnextflowpipeline)"; echo "-m Path to MultiQC_SSDS_Rev1 directory (default : /home/${USER}/work/ssdsnextflowpipeline/accessoryFiles/SSDS/MultiQC_SSDS_Rev1)"; echo "-e Name for Nextflow conda environment (default : nextflow_dev)"; echo "-c Name for multiqc0.7-dev conda environment (default : multiqc_dev)"; exit 0;;
        p) PIPELINE_DIRECTORY=${OPTARG};if [ ! -d ${PIPELINE_DIRECTORY} ]; then echo "Directory ${PIPELINE_DIRECTORY} not found!" ; exit 0; fi;;
        m) MULTIQC_DEV_DIRECTORY=${OPTARG};if [ ! -d ${MULTIQC_DEV_DIRECTORY=} ]; then echo "Directory ${MULTIQC_DEV_DIRECTORY} not found!" ; exit 0; fi;;
        e) MAIN_ENV_NAME=${OPTARG};;
	c) MULTIQC_DEV_NAME=${OPTARG};;
	f) FLAG=${OPTARG}; if [ $FLAG != 3 ] &&[ $FLAG != 2 ] && [ $FLAG != 1 ] && [ $FLAG != 0 ] ; then echo "-f value should be 0 ; 1 ; 2 or 3 see help"; exit 0; fi;;
    esac
done

#Go to pipeline directory
cd ${PIPELINE_DIRECTORY}

eval "$(conda shell.bash hook)"

if [ $FLAG == 1 ] || [ $FLAG == 2 ] ;
then
## Create conda environment to install nextflow
conda create -y -n ${MAIN_ENV_NAME}
conda activate ${MAIN_ENV_NAME}
conda install -c bioconda nextflow=20.04.1
conda deactivate 

elif [ $FLAG == 3 ] || [ $FLAG == 2 ] ;
then
## Create conda environment for multiQC0.7 dev library
#conda create -y -p ${PIPELINE_DIRECTORY}/${MULTIQC_DEV_NAME} 
#conda activate ${PIPELINE_DIRECTORY}/${MULTIQC_DEV_NAME}
#conda env config vars set PYTHONPATH="${MULTIQC_DEV_DIRECTORY}/lib/python2.7/site-packages/"
#conda activate ${PIPELINE_DIRECTORY}/${MULTIQC_DEV_NAME}

## Install Multiqc_dev python library in the dedicated conda env
#cd ${MULTIQC_DEV_DIRECTORY}
#python setup.py build
#python setup.py install --prefix "${MULTIQC_DEV_DIRECTORY}"
#echo ${PIPELINE_DIRECTORY}
#echo ${MULTIQC_DEV_DIRECTORY}
#conda deactivate
echo "Multiqc dev is not possible at the moment, sorry."
fi ;

