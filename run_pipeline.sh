#!/bin/bash

#Conda parameters
CONDA_ENV="SSDSnextflowPipeline"

#Activate conda environement
eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV

#Pipeline parameters
WORKING_DIRECTORY=/home/demassyie/work
PROJECT_NAME="runtest_SSDSnextflowPipeline"
GENOME="mm10"
R1LEN="36"
R2LEN="40"
FQ1=$NXF_PIPEDIR/tests/fastq/ssdsLong.100k.R1.fastq
FQ2=$NXF_PIPEDIR/tests/fastq/ssdsLong.100k.R2.fastq

#Create working directory if not existing
mkdir -p ${WORKING_DIRECTORY}/${PROJECT_NAME}
cd ${WORKING_DIRECTORY}/${PROJECT_NAME}

#Run nextflow
nextflow run -c ${NXF_PIPEDIR}/conf/shenron.config \
    ${NXF_PIPEDIR}/SSDSPipeline_1.8.nf \
    --fq1 ${FQ1} \
    --fq2 ${FQ2} \
    --r1Len ${R1LEN} \
    --r2Len ${R2LEN} \
    --genome ${GENOME} \
    --name ${PROJECT_NAME} \
    --outdir ${WORKING_DIRECTORY}/${PROJECT_NAME} \
    -resume

#deactivate conda environment
conda deactivate
