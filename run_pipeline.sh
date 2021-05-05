#!/bin/bash

#Run ssds + callpeaks nextflow pipeline

#Set slurm environment variables
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

#Singularity settings
#Set environmental variable to mount the non canonical directories on IGH cluster
export SINGULARITY_BINDPATH="/work,/poolzfs"

#Pipeline default parameters
ANALYSIS_NAME="SSDS_pipeline"
PIPELINE_DIRECTORY="/home/${USER}/work/ssdsnextflowpipeline"
OUTPUT_DIRECTORY="/home/${USER}/work/results/${ANALYSIS_NAME}.outdir"
CONF="${PIPELINE_DIRECTORY}/conf/igh.config"
CENV="/home/${USER}/work/bin/miniconda3/envs/nextflow_dev"
JOBNAME='SSDS_main'
INPUT=""
OPTIONS="-profile conda "
TEST="0"

#Get command line arguments
while getopts hp:o:c:a:i:n:y:t: flag
do
	case "${flag}" in
		h) echo ""; echo "Usage: bash `basename $0` -i input_file [options] "; echo "Options : "; echo "-h display help message"; echo "-i Absolute path to input csv file (REQUIRED except in case of run test on test dataset) "; echo "-p Absolute path to ssds nextflow pipeline base directory (default : ${PIPELINE_DIRECTORY})"; echo "-o Absolute path to output directory (default : ${OUTPUT_DIRECTORY})"; echo "-c Absolute path to IGH cluster configuration file (default : ${CONF})"; echo "-a Absolute path to conda environment for nextflow (default : ${CENV})"; echo "-n Analysis name (default : ${ANALYSIS_NAME})"; echo "-y Optional arguments for the pipeline (for example \"--with_control --no_multimap --trim_cropR1 50 --trim_cropR2 50\" ;  default : \"${OPTIONS}\")"; echo "-t set to 1 if running pipeline on test data located in ${PIPELINE_DIRECTORY}/tests/fastq (default : ${TEST})"; echo ""; exit 0;;
		p) PIPELINE_DIRECTORY=${OPTARG};if [ ! -d ${PIPELINE_DIRECTORY} ]; then echo "Directory ${PIPELINE_DIRECTORY} not found!" ; exit 0; fi;;
		o) OUTPUT_DIRECTORY=${OPTARG};;
		c) CONF=${OPTARG};if [ ! -f ${CONF} ]; then echo "File ${CONF} not found!" ; exit 0; fi;;
		a) CENV=${OPTARG};if [ ! -d ${CENV} ]; then echo "Environment ${CENV} not found!" ; exit 0; fi;;
		i) INPUT=${OPTARG};if [ ! -f ${INPUT} ]; then echo "File ${INPUT} not found!" ; exit 0; fi;;
		y) OPTIONS=${OPTARG};;
		n) ANALYSIS_NAME=${OPTARG};;
		t) TEST=${OPTARG};;
		\? ) echo "Unknown option: -$OPTARG" >&2; exit 1;;
        	:  ) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
        	*  ) echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
	esac
done

#Option -i is mandatory
if [[ "x" == "x$INPUT" && ${TEST} == "0" ]] ; then
  echo "-i [input csv file] is required"
  exit
fi

#If running pipeline on test dataset, create relevant input csv file
if [ ${TEST} == "1" ]; then
	echo "Running SSDS pipeline on test dataset. Creating input file..."
	echo "group,replicate,fastq_1,fastq_2,antibody,control" > ${PIPELINE_DIRECTORY}/tests/fastq/input.csv
	echo "TEST_IP,1,${PIPELINE_DIRECTORY}/tests/fastq/ssdsLong.100k.R1.fastq,${PIPELINE_DIRECTORY}/tests/fastq/ssdsLong.100k.R2.fastq,antiDMC1," >> ${PIPELINE_DIRECTORY}/tests/fastq/input.csv
	INPUT=${PIPELINE_DIRECTORY}/tests/fastq/input.csv
	if [ -f ${INPUT} ]; then 
		echo "Input file created ; check ${PIPELINE_DIRECTORY}/tests/fastq/input.csv."
	else
		echo "Input file for test data cannot be created."
		exit 0
	fi
fi 

#Create output directory if not existing
mkdir -p ${OUTPUT_DIRECTORY}
cd ${OUTPUT_DIRECTORY}

#Activate conda environment
eval "$(conda shell.bash hook)"
conda activate ${CENV}

#Run pipeline
echo "Running SSDS pipeline from ${PIPELINE_DIRECTORY} on ${INPUT##*/} data within ${CENV##*/} conda environment. Check output directory ${OUTPUT_DIRECTORY}"

#Cheat lines for dev, do not use
#OPTIONBASE="--with_control --satcurve false --kbrick_bigwig false -with-tower --bigwig_profile T12rep  --genome mm10 --no_multimap --trim_cropR1 50 --trim_cropR2 50 --nb_replicates 2 --with_ssds_multiqc --multiqc_dev_conda_env /work/demassyie/bin/miniconda2/envs/SSDSnextflowPipeline -resume"
#OPTIONBASE="--genome mm10 --no_multimap --trim_cropR1 50 --trim_cropR2 50 --binsize 25  -resume"
OPTIONBASE=""

sbatch -p computepart -J ${JOBNAME} --export=ALL -n 1 --mem 7G -t 5-0:0  \
--wrap "export MKL_NUM_THREADS=1 ; export NUMEXPR_NUM_THREADS=1 ; export OMP_NUM_THREADS=1 ; \
nextflow run ${PIPELINE_DIRECTORY}/main.nf -c ${CONF} --name ${ANALYSIS_NAME} --outdir ${OUTPUT_DIRECTORY} --inputcsv ${INPUT} ${OPTIONBASE} ${OPTIONS}"

#Deactivate conda environment
conda deactivate

