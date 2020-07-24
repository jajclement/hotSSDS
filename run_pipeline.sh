#!/bin/bash

#Conda parameters
CONDA_ENV="SSDSnextflowPipeline"

#Activate conda environement
eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV

#Pipeline parameters
WORKING_DIRECTORY=/home/demassyie/work
#PROJECT_NAME="runtest_SSDSnextflowPipeline"
PROJECT_NAME="results_SSDSnextflowPipeline_17144FL-05-01-02_S2"
GENOME="mm10"
R1LEN="36"
R2LEN="40"
#FQ1=${NXF_PIPEDIR}/tests/fastq/ssdsLong.100k.R1.fastq
#FQ2=${NXF_PIPEDIR}/tests/fastq/ssdsLong.100k.R2.fastq
#FQ1=${NXF_PIPEDIR}/input_data/raw_data/17144FL-05-01-01_S1_L001_R1_001.fastq.gz
#FQ2=${NXF_PIPEDIR}/input_data/raw_data/17144FL-05-01-01_S1_L001_R2_001.fastq.gz
FQ1=${NXF_PIPEDIR}/input_data/raw_data/17144FL-05-01-02_S2_L001_R1_001.fastq.gz
FQ2=${NXF_PIPEDIR}/input_data/raw_data/17144FL-05-01-02_S2_L001_R2_001.fastq.gz
DAG=${WORKING_DIRECTORY}/${PROJECT_NAME}/dag.png

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
    -resume \
    -with-dag ${DAG} \
    -with-tower

#deactivate conda environment
conda deactivate
