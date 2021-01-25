#!/bin/bash

#Run ssds + callpeaks nextflow pipeline

#Set slurm environment variables
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

#Default parameters
ANALYSIS_NAME="SSDS_pipeline_libsize"
PIPELINE_DIRECTORY="/home/${USER}/work/ssdsnextflowpipeline"
OUTPUT_DIRECTORY="/home/${USER}/work/results/${ANALYSIS_NAME}.outdir"
CONF="${PIPELINE_DIRECTORY}/conf/igh.config"
CENV=/home/${USER}/work/bin/miniconda3/envs/nextflow-dev
JOBNAME='SSDS_main'
INPUT=""
OPTIONS="-profile conda " # -with-tower"

#Get command line arguments
while getopts hp:o:c:a:i:n:y: flag
do
	case "${flag}" in
		h) echo ""; echo "Usage: bash `basename $0` -i input_file [options] "; echo "Options : "; echo "-h display help message"; echo "-i REQUIRED Absolute path to input csv file"; echo "-p Path to ssds nextflow pipeline  base directory (default : ${PIPELINE_DIRECTORY})"; echo "-o Path to output directory (default : ${OUTPUT_DIRECTORY})"; echo "-c Path to IGH cluster configuration file (default : ${CONF})"; echo "-a Path to conda environment for nextflow (default : ${CENV})"; echo "-n Analysis name (default : ${ANALYSIS_NAME})"; echo "-y Optional arguments for the pipeline (for example \"--with_control --no_multimap --trim_cropR1 50 --trim_cropR2 50\" ;  default : \"${OPTIONS}\")"; echo ""; exit 0;;
		p) PIPELINE_DIRECTORY=${OPTARG};if [ ! -d ${PIPELINE_DIRECTORY} ]; then echo "Directory ${PIPELINE_DIRECTORY} not found!" ; exit 0; fi;;
		o) OUTPUT_DIRECTORY=${OPTARG};;
		c) CONF=${OPTARG};if [ ! -f ${CONF} ]; then echo "File ${CONF} not found!" ; exit 0; fi;;
		a) CENV=${OPTARG};if [ ! -d ${CENV} ]; then echo "Environment ${CENV} not found!" ; exit 0; fi;;
		i) INPUT=${OPTARG};if [ ! -f ${INPUT} ]; then echo "File ${INPUT} not found!" ; exit 0; fi;;
		y) OPTIONS=${OPTARG};;
		n) ANALYSIS_NAME=${OPTARG};;
		\? ) echo "Unknown option: -$OPTARG" >&2; exit 1;;
        	:  ) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
        	*  ) echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
	esac
done

#Option -i is mandatory
if [ "x" == "x$INPUT" ]; then
  echo "-i [input csv file] is required"
  exit
fi

#Create output directory if not existing
OUTPUT_DIRECTORY="/home/${USER}/work/results/${ANALYSIS_NAME}.outdir"
mkdir -p ${OUTPUT_DIRECTORY}
cd ${OUTPUT_DIRECTORY}

#Activate conda environment
eval "$(conda shell.bash hook)"
conda activate ${CENV}

#Run pipeline
echo "Running pipeline nf-core chipseq ${PIPELINE_DIRECTORY##*/} on ${INPUT##*/} data within ${CENV##*/} conda environment. Check output directory ${OUTPUT_DIRECTORY}"

#Cheat lines for dev
OPTIONBASE="--nb_replicates 2 --bigwig_profile minimal2 --satcurve false --genome mm10 --no_multimap  --trim_cropR1 50 --trim_cropR2 50 --binsize 25  --with_ssds_multiqc --multiqc_dev_conda_env /work/demassyie/bin/miniconda2/envs/SSDSnextflowPipeline -resume"
#OPTIONBASE="--genome mm10 --no_multimap --trim_cropR1 50 --trim_cropR2 50 --binsize 25  -resume"
#OPTIONBASE=""

sbatch -p computepart -J ${JOBNAME} --export=ALL -n 1 --mem 7G -t 5-0:0 --mem-per-cpu=1000 \
--wrap "export MKL_NUM_THREADS=1 ; export NUMEXPR_NUM_THREADS=1 ; export OMP_NUM_THREADS=1 ; nextflow run ${PIPELINE_DIRECTORY}/main.nf -c ${CONF} --name ${ANALYSIS_NAME} --outdir ${OUTPUT_DIRECTORY} --inputcsv ${INPUT} ${OPTIONBASE} ${OPTIONS}"

#Deactivate conda environment
conda deactivate

