#!/bin/bash

WORKING_DIRECTORY=/home/demassyie/work/SSDSnextflowPipeline
GENOMES_DIR_LOCAL=${WORKING_DIRECTORY}/input_data/genomes/mm10
INDEX_DIR_LOCAL=${GENOMES_DIR_LOCAL}/BWAIndex/version0.7.10
GENOMES_DIR_CLUSTER=/poolzfs/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta
INDEX_DIR_CLUSTER=/poolzfs/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta
DATA_DIR=/home/demassyie/work/raw_data/Nore_may2019/17144-05
DATA_DIR_PROJECT=${WORKING_DIRECTORY}/input_data/raw_data

mkdir -p ${WORKING_DIRECTORY}
cd ${WORKING_DIRECTORY}

echo "working dir : ${WORKING_DIRECTORY}"
echo "genome dir cluster : ${GENOMES_DIR_CLUSTER}"
echo "genome dir local : ${GENOMES_DIR_LOCAL}"
echo "index dir local : ${INDEX_DIR_LOCAL}"
echo "index dir cluster : ${INDEX_DIR_CLUSTER}"
echo "data dir : ${DATA_DIR}"
echo "data dir project : ${DATA_DIR_PROJECT}"

mkdir -p ${GENOMES_DIR_LOCAL}
mkdir -p ${INDEX_DIR_LOCAL}
mkdir -p ${DATA_DIR_PROJECT}

#ln -s ${GENOMES_DIR_CLUSTER}/genome.fa ${GENOMES_DIR_LOCAL}/genome.fa
#ln -s ${GENOMES_DIR_CLUSTER}/genome.fa.fai ${GENOMES_DIR_LOCAL}/genome.fa.fai
#ln -s ${GENOMES_DIR_CLUSTER}/genome.dict ${GENOMES_DIR_LOCAL}/genome.dict
#ln -s ${INDEX_DIR_CLUSTER}/* ${INDEX_DIR_LOCAL}/
#ln -s ${DATA_DIR}/* ${DATA_DIR_PROJECT}/

