#!/bin/bash

#Run ssds + callpeaks nextflow pipeline

#Singularity settings
#Set environmental variable to mount the non canonical directories on IGH cluster
export SINGULARITY_BINDPATH=""

#Pipeline default parameters
ANALYSIS_NAME="SSDS_pipeline"
PIPELINE_DIRECTORY="${DATAWORK}/ssdsnextflowpipeline"
BASE_DIRECTORY="${DATAWORK}/results"
CONF="${PIPELINE_DIRECTORY}/conf/cluster.config"
GENOME_PROFILE="${PIPELINE_DIRECTORY}/conf/mm10.json"
now=`date +"%FT%H%M%S"`
INPUT=""
OPTIONS=""
TEST="0"
FORCE="0"
PARAMS_FILE="${PIPELINE_DIRECTORY}/conf/mm10.config"
TOWER_TOKEN="none"

#Get command line arguments
while getopts hp:b:n:c:i:o:w:t:f:g: flag
do
	case "${flag}" in
		h) echo ""; echo "Usage: bash `basename $0` -i input_file [options] "; \
		   echo "Options : "; echo "-h display help message"; \
		   echo "-i Absolute path to input csv file (REQUIRED except in case of run test on test dataset) "; \
		   echo "-g Absolute path to the genome config file (default : ${GENOME_PROFILE})" ; \
		   echo "-p Absolute path to ssds nextflow pipeline base directory (default : ${PIPELINE_DIRECTORY})"; \
		   echo "-b Absolute path to base directory where the output directory will be created (default : ${BASE_DIRECTORY})"; \
		   echo "-n Analysis name (default : ${ANALYSIS_NAME}) INFO : by default, this parameter will match the --name option in nextflow command line"; \
		   echo "-c Absolute path to IGH cluster configuration file (default : ${CONF})"; \
		   echo "-o Optional arguments for the pipeline (for example \"--with_control --no_multimap --trim_cropR1 50 --trim_cropR2 50\" ;  default : \"${OPTIONS}\")"; \
		   echo "-w Valid Nextflow Tower token (default : ${TOWER_TOKEN} ; if not none, then the option -with-tower has to be added in -o parameter))"; \
		   echo "-t set to 1 if running pipeline on test data located in ${PIPELINE_DIRECTORY}/tests/fastq (default : ${TEST})"; \
		   echo "-f set to 1 to force pipeline to run without checking resume/output directory (default : ${FORCE})" ; \
		   echo "INFO : the output directory will be located in the base directory and will be named after the analysis name parameter with the .outdir suffix (default ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir)"; \
		   echo ""; exit 0;;
		p) PIPELINE_DIRECTORY=${OPTARG};if [ ! -d ${PIPELINE_DIRECTORY} ]; then echo "Directory ${PIPELINE_DIRECTORY} not found!" ; exit 0; fi;;
		b) BASE_DIRECTORY=${OPTARG};;
		n) ANALYSIS_NAME=${OPTARG};;
		g) GENOME_PROFILE=${OPTARG};if [ ! -f ${GENOME_PROFILE} ]; then echo "File ${GENOME_PROFILE} not found!" ; exit 0; fi;;
		c) CONF=${OPTARG};if [ ! -f ${CONF} ]; then echo "File ${CONF} not found!" ; exit 0; fi;;
		i) INPUT=${OPTARG};if [ ! -f ${INPUT} ]; then echo "File ${INPUT} not found!" ; exit 0; fi;;
		o) OPTIONS=${OPTARG};;
		w) TOWER_TOKEN=${OPTARG};;
		t) TEST=${OPTARG};;
		f) FORCE=${OPTARG};;
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

#Create base directory if not existing
echo "Checking base directory..."
mkdir -p ${BASE_DIRECTORY}

#Test if a project with the same name already exists and if the option -resume is properly set
if [[ $FORCE == 0 ]] 
then
	echo "Checking output directory."
	if [[ -d ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir && $OPTIONS == *"-resume"* ]]
	then
		echo "WARNING : Output directory already exists for project ${ANALYSIS_NAME} in ${BASE_DIRECTORY}, are you sure you want to resume this run ? (y/n)"
		read answer1
		case $answer1 in
			[yYoO]*) echo "Ok; run will be resumed.";;
			[nN]*) echo "Ok; quitting. Bye, see you soon !"; exit 0;;
			*) echo "ABORT. Please enter y (yes) or n (no) next time. Bye, see you soon !"; exit 1;;
		esac
	elif [[ -d ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir && ! $OPTIONS == *"-resume"* ]]
	then
		echo "WARNING : Output directory already exists for project ${ANALYSIS_NAME} in ${BASE_DIRECTORY}, are you sure you want to run the pipeline from the beginning (this will erase previous results from this project) ? (y/n)"
		read answer2
		case $answer2 in
			[yYoO]*) echo "Ok; previous results located in ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir will be erased.";;
			[nN]*) echo "Ok; quitting. Consider using another analysis name (-n argument) or using the option -resume in the -o argument next time. Bye, see you soon !"; exit 0;;
			*) echo "ABORT. Please enter y (yes) or n (no) next time. Bye, see you soon !"; exit 1;;
		esac
	else
		echo "Ok, creating ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir output project directory."
		mkdir ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir
	fi
else
	echo "Forcing the creation of ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir output project directory."
	mkdir -p ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir
fi

#Go to the output directory
cd ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir

#Create output directory for log files
mkdir -p ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir/slurm

#Check pipeline directory
echo "Checking pipeline directory..."
if [ ! -d ${PIPELINE_DIRECTORY} ];
then echo "ABORT : pipeline directory ${PIPELINE_DIRECTORY} not found! Check -p argument. Bye, see you soon !" ; exit 0;
else echo "Ok.";
fi

#Checking configuration files
echo "Checking configuration files..."
if [ ! -f ${CONF} ] ;
then echo "ABORT : configuration file ${CONF} not found! Check -c argument. Bye, see you soon !" ; exit 0;
else echo "Ok.";
fi
if [ ! -f ${GENOME_PROFILE} ];
then echo "ABORT : configuration file ${GENOME_PROFILE} not found! Check -g argument. Bye, see you soon !" ; exit 0; 
else echo "Ok.";
fi 

#Run pipeline
echo "Run pipeline !"
if [ ${TEST} == "0" ]; then
echo "Running SSDS pipeline from ${PIPELINE_DIRECTORY} on ${INPUT##*/} data. Check output directory ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir/"
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

#Cheat lines for dev, do not use
#OPTIONBASE="--with_control --satcurve false --kbrick_bigwig false -with-tower --bigwig_profile T12rep  --genome mm10 --no_multimap --trim_cropR1 50 --trim_cropR2 50 --nb_replicates 2 --with_ssds_multiqc --multiqc_dev_conda_env /work/demassyie/bin/miniconda2/envs/SSDSnextflowPipeline -resume"
#OPTIONBASE="--genome mm10 --no_multimap --trim_cropR1 50 --trim_cropR2 50 --binsize 25  -resume"
OPTIONBASE=""

CMD="nextflow run ${PIPELINE_DIRECTORY}/main.nf -c ${CONF} -params-file ${GENOME_PROFILE} --name ${ANALYSIS_NAME} --outdir ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir --inputcsv ${INPUT} ${OPTIONBASE} ${OPTIONS}"

# Launch PBS run script to run the pipeline on HPC cluster using pbs pro
echo
qsub -v CMDLINE="${CMD}",TOWER_TOKEN=${TOWER_TOKEN},OUTDIR=${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir ${PIPELINE_DIRECTORY}/run_main.pbs
echo

