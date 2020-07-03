#!/bin/bash

CONDA_ENV="SSDSnextflowPipeline"


#Activate conda environement
eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV

#Run nextflow
nextflow run -c $NXF_PIPEDIR/conf/shenron.config \
    $NXF_PIPEDIR/SSDSPipeline_1.8.nf \
#    --fq1 $NXF_PIPEDIR/tests/fastq/ssdsLong.100k.R1.fastq \
#    --fq2 $NXF_PIPEDIR/tests/fastq/ssdsLong.100k.R2.fastq \
#    --r1Len 36 \
#    --r2Len 40 \
#    --genome mm10 \
#    --name testSSDS1.8 \
#    --outdir SSDS1.8_test

#deactivate conda environment
conda deactivate
