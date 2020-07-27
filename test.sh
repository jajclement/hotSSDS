#!/bin/bash

WORKDIR='/home/demassyie/work/SSDSnextflowPipeline'
SCRIPT=${WORKDIR}/test.nf
CONF=${WORKDIR}/test.conf
CONDA_ENV='SSDSnextflowPipeline'

eval "$(conda shell.bash hook)"
conda activate SSDSnextflowPipeline

nextflow run ${SCRIPT} -c ${CONF} #-resume 

conda deactivate
