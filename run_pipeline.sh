#!/bin/bash

WORKDIR='/home/demassyie/work/SSDSnextflowPipeline'
SCRIPT=${WORKDIR}/main.nf
CONF=${WORKDIR}/conf/igh.config
CONDA_ENV='SSDSnextflowPipeline'

eval "$(conda shell.bash hook)"
#conda activate nextflow-dev
conda activate ${CONDA_ENV}

nextflow run ${SCRIPT} -c ${CONF} \
	--fqdir "${WORKDIR}/tests/fastq/*{R1,R2}.fastq" \
	--name "debug" \
	-resume
	#--fqdir "${WORKDIR}/input_data/raw_data/*{1,2}_001.fastq.gz" \
	#--name "17144FL_dataset_crop50" -resume \
	#-profile conda \
	#-resume

conda deactivate
