#!/bin/bash

#ANALYSIS_NAME="17144FL_dataset_crop50"
#ANALYSIS_NAME="DATA_KHIL_2012"
ANALYSIS_NAME="TEST"
WORKDIR="/home/${USER}/work"
#INPUT="${WORKDIR}/SSDSnextflowPipeline/input_data/raw_data/*{R1,R2}_001.fastq.gz"
#INPUT=['SRR391617','SRR391618','SRR391619','SRR391620','SRR391621','SRR391622','SRR391623','SRR391624','SRR391625','SRR391626','SRR391627','SRR391628','SRR391629','SRR391630','SRR391631','SRR391632','SRR391633','SRR391634','SRR391635','SRR391636','SRR391637','SRR391638','SRR391639','SRR391640']
#INPUT='SRR391633'
INPUT="${WORKDIR}/SSDSnextflowPipeline/tests/fastq/*{R1,R2}.fastq"

SCRIPT="${WORKDIR}/SSDSnextflowPipeline/main.nf"
CONF="${WORKDIR}/SSDSnextflowPipeline/conf/igh.config"
CONDA_ENV='nextflow-dev'
JOBNAME='SSDS_main'

eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}

sbatch -p computepart -J ${JOBNAME} --export=ALL --mem 5G -t 5-0:0 --mem-per-cpu=1000 \
	--wrap "nextflow run ${SCRIPT} -c ${CONF} --fqdir ${INPUT} --name ${ANALYSIS_NAME} --trim_cropR1 36 --trim_cropR1 40 -profile conda -resume"

conda deactivate

