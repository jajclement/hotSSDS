#!/bin/bash

ANALYSIS_NAME="TEST"
INPUT="tests/fastq/*{R1,R2}.fastq"
GENOME="mm10"
JOBNAME='SSDS_main'

eval "$(conda shell.bash hook)"
conda activate nextflow-dev

sbatch -p computepart -J ${JOBNAME} --export=ALL --mem 5G -t 5-0:0 --mem-per-cpu=1000 \
	--wrap "nextflow run mains.nf -c conf/igh.config --fqdir ${INPUT} --name ${ANALYSIS_NAME} --genome ${GENOME}"

conda deactivate

