## **SSDS nextflow pipeline version 2.0**
### **Parse, Align and Call peaks from SSDS data**
**\/!\ Work in progress /!\ **


### **Welcome**
**This pipeline is based on  [SSDS pipeline](https://github.com/kevbrick/SSDSnextflowPipeline) and [SSDS call peaks pipeline](https://github.com/kevbrick/callSSDSpeaks) by Kevin Brick, updated and adapted to IGH cluster.**
See [initial paper](https://genome.cshlp.org/content/22/5/957.long) and [technical paper](https://www.sciencedirect.com/science/article/pii/S0076687917303750?via%3Dihub).
The pipeline uses [Nextflow]( https://www.nextflow.io/) > 20.04.1
Briefly, the update from SSDS pipeline version 1.8_NF included **conda profile**, input modification, callpeaks and IDR procedure addition, and global nextflow homogeneisation.

*Please report any bug occured during the installation to pauline.auffret at igh dot cnrs dot fr*

The pipeline is composed of 25 processes :
*PROCESS 1 : check_design (CHECK INPUT DESIGN FILE)
*PROCESS 2 : makeScreenConfigFile (MAKE CONFIGURATION FILE FOR FASTQSCREEN)
*PROCESS 3 : trimming (USE TRIMMOMATIC OR TRIM-GALORE TO QUALITY TRIM, REMOVE ADAPTERS AND HARD TRIM SEQUENCES)
*PROCESS 4 : fastqc (QUALITY CONTROL ON RAW READS USING FASTQC AND FASTQSCREEN)
*PROCESS 5 : bwaAlign (USE BWA AND CUSTOM BWA (BWA Right Align) TO ALIGN SSDS DATA)
*PROCESS 6 : filterBam (MARK DUPLICATES, REMOVE SUPPLEMENTARY ALIGNMENTS, SORT AND INDEX)
*PROCESS 7 : parseITRs (PARSE BAM FILES TO GET THE 5 DIFFERENT SSDS TYPES : SSTYPE1, SSTYPE2, DS, DSSTRICT, UNCLASSIFIED)
*PROCESS 8 : makeBigwig (GENERATE NORMALIZED BIGWIG FILES FOR T1 AND T1+T2 BED FILES)
*PROCESS 9 : makeBigwigReplicates (optional, GENERATE NORMALIZED BIGWIG FILES FOR MERGED REPLICATES T1 AND T1+T2 BED FILES)
*PROCESS 10 : makeDeeptoolsBigWig (optional, GENERATES BIGWIG FILES, COVERAGE AND CUMULATIVE COVERAGE PLOTS)
*PROCESS 11 : toFRBigWig (optional, GENERATES FWD/REV BIGWIG FILES)
*PROCESS 12 : samStats (GENERATES SAMSTATS REPORTS)
*PROCESS 13 : makeSSreport (COMPUTE STATS FROM SSDS PARSING)
*PROCESS 14 : ssds_multiqc (MAKE MULTIQC REPORT FOR SSDS FILES)
*PROCESS 15 : makeFingerPrint (MAKE DEEPTOOLS FINGERPRINT PLOTS)
*PROCESS 16 : shufBEDs (T1 BED SHUFFLING)
*PROCESS 17 : callPeaks (PEAK CALLING ON T1 WITH MACS2)
*PROCESS 18 : createPseudoReplicates (optional, CREATES ALL PSEUDOREPLICATES AND POOL FOR IDR ANALYSIS)
*PROCESS 19 : callPeaksForIDR (optional, CALL PEAKS WITH MAC2 ON ALL REPLICATES AND PSEUDO REPLICATES)
*PROCESS 20 : IDRanalysis (optional, PERFORM IDR ANALYSIS ON 4 PAIRS OF REPLICATES OR PSEUDOREPLICATES
*PROCESS 21 : IDRpostprocess (optional, IDR PEAKS POST PROCESSING)
*PROCESS 22 : normalizePeaks (CENTER AND NORMALIZE PEAKS)
*PROCESS 23 ; mergeFinalPeaks (MERGE PEAKS FROM REPLICATES)
*PROCESS 24 : makeSatCurve (optional, CREATE SATURATION CURVE)
*PROCESS 25 : general_multiqc (GENERATES GENERAL MULTIQC REPORT)

## **How to run the pipeline on IGH cluster**
### Requirements
Get [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) if not installed on your system. First, download the installation script (for example in /home/${USER}/work/bin directory) :
````sh
mkdir -p /home/${USER}/work/bin
cd /home/${USER}/work/bin
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
````
Then execute the installation script and follow the prompts on the installer screens :
````sh
bash Miniconda3-latest-Linux-x86_64.sh
````

### 1. Get the pipeline and set up conda environment
First you need to download the pipeline in your working directory (in the following instructions /home/${USER}/work/ will refer to your working directory. Please substitute with the according path if different) :
````sh
cd /home/${USER}/work
git clone https://gitlab.igh.cnrs.fr/pauline.auffret/ssdsnextflowpipeline.git
cd ssdsnextflowpipeline
````
Then install pipeline conda environment through sbatch script : 
```` 
sbatch -p computepart -J "install_conda_env" --export=ALL --mem 5G -t 5-0:0 --mem-per-cpu=1000 --wrap "bash src/install_pipeline.sh"
```` 
This will create 1 conda environment : **nextflow-dev**.
Please use ``bash src/install_pipeline.sh -h`` to see details.

### 2. Pipeline configuration 
There are currently 2 configuration files :
- ````./conf/igh.config```` contains cluster resources requirements & reference genomes info. You don't need to edit this file unless you want to custom the requirements for CPU/memory usage and compute queue (see [IGH cluster documentation](https://kojiki.igh.cnrs.fr/doku.php?id=cluster,))
- ````./nextflow.config```` contains default pipeline parameters. You *can* edit this file but it is recommended to use parameters in nextflow command line through ````run_pipeline.sh```` script to use your own parameters (this will overwrite the default configuration).
If you are running the pipeline on an other computing cluster, you need to specify the relevant configuration file.


### 3. Input data
The pipeline will only process **paired-end data**.
The input data must be described in an input csv file with the 6 following fields (see template in ``/home/${USER}/work/ssdsnextflowpipeline//tests/fastq/input.csv``)
````
group,replicate,fastq_1,fastq_2,antibody,control
DMC1-chip,1,/work/${USER}/data/SRR1035576_R1.fastq.gz,/work/${USER}/data/SRR1035576_R2.fastq.gz,antiDMC1,Input
DMC1-chip,2,/work/${USER}/data/SRR1035577_R1.fastq.gz,/work/${USER}/data/SRR1035577_R2.fastq.gz,antiDMC1,Input
Input,1,/work/${USER}/data/SRR1035578_R1.fastq.gz,/work/${USER}/data/SRR1035578_R2.fastq.gz,,,
````

The reference genome should be in the ``/poolzvs/genomes`` directory on IGH cluster. Currently, the available genomes are mm10, hg19, hg38, sacCer2, sacCer3, dm3, dm6.

You can use ````-- genome 'mm10'```` and you won't need to worry about the other genome parameters like fasta path etc.

But if you want to use another reference, you will need to set : 
- path to genome ````--genomedir /path/to/genome````
- path to fasta file. Indexes for BWA SHOULD EXIST in the same directory ````--genome_fasta /path/to/genome.fa````
- the name of the genome ````--genome_name mm11````
- the path to the fai index ````--fai /path/to/genome.fai.fai````

### 4. Run the pipeline !
To see all the available parameters, please run :
````
cd /home/${USER}/work/ssdsnextflowpipeline
conda activate nextflow-dev
nextflow run main.nf --help
````

The main parameters you need to set are :
* ``--inputcsv`` : path to the input csv file, see section 3
* ``--profile conda`` : Run within conda environment (only available option at the moment)
* ``--genome`` : the reference genome, see section 3
* ``--name`` : analysis name, e.g. "SSDS_SRA5678_DMC1"
* ``--outdir`` : path to output directory
* ``--with_control`` (true/false) : use input control files, see section 3
* ``--no_multimap`` (true/false) : remove multimappers from bam files
* ``--nb_replicates`` : number of biological replicates you are running with (maximum 2)
* ``--bigwig_profile`` : indicates which bigwig to generate (T1 ; T12 ; T1rep or T12rep)
* ``--with_idr`` (true/false) : run IDR analysis (if nb_replicates=2)
* ``--satcurve`` (true/false) : plot saturation curve



You can either run the pipeline directly through the command line :
````
cd /home/${USER}/work/ssdsnextflowpipeline
conda activate nextflow-dev
sbatch -p computepart -J SSDSnextflowPipeline --export=ALL --mem 5G -t 5-0:0 --mem-per-cpu=1000 --wrap "nextflow run main.nf -c conf/igh.config --inputcsv /path/to/the/input/csv/file --name your_analysis_name --genome mm10 --profile conda --with_control --nb_replicates 2"
````
or use ``run_pipeline.sh`` script :
````
bash run_pipeline.sh -i inputfile.csv -y "--name your_analysis_name --genome mm10 --profile conda --with_control --nb_replicates 2"
````
See the available options before use with 
````
bash run_pipeline.sh -h
````

**You can use ``-resume`` option (Nextflow native options) to prevent the entire workflow to be relaunch in case you need to relaunch an aborted workflow** 




**Nextflow Tower**
You can use ``-with-tower`` option to monitor your jobs through [nextflow tower web interface](https://tower.nf/). 
You first need to sign in to get your key, then add it to your parameters with the ``--tower-token 'yourkey'`` option.

### 5. Output
The main output folder is specified through the ````--outdir```` parameter.
This folder will contain the following directories :
* **pipeline_info** with input csv files rearranged for the pipeline
* **trim_fastqc** with fastQC reports for trimmed fastq files
* **trim_fastq** with trimmed fastq files
* **fastqscreen** with plots from fastqscreen 
* **conf.fqscreen** is the input file for fastqscreen
* **raw_fastqc** with fastQC reports for raw fastq files
* **bam** with bam and bai files before filtering
* **filterbam** with filtered bam and bai files 
* **parse_itr** with bam, bai and bed files parsed from filtered bam files to 5 subtypes : ssT1, ssT2, ds, ds_strict, unclassified
* **samstats** with samtools statistics reports for bam files
* **bed_shuffle** with shuffled T1 bed files before callpeaks
* **multiqc** with multiqc quality reports
* **bigwig** with bigwig and bedgraph files
* **pseudo_replicates** with pseudoreplicates bed files if IDR analysis is run
* **peaks** with peaks bed files from callpeaks
* **saturation_curve** with plots and bed peak files for saturation curve of samples
* **idpeaks** with peaks bed files after IDR if IDR analysis is run
* **idrresults** with IDR final peak bed files
* **idrQC** with reports from IDR analysis
* **finalpeaks** with final peaks set from the pipeline, to use for downstream analysis
* **normpeaks** with recentered and normalized bedgraphs files after callpeaks
* **work** is the working directory for Nextflow
* **fingerprint** with fingerprint plot for filtered bam files
* **nxfReports** with nextflow execution reports and diagram


The execution reports are in ``nxfReports`` folder created in your output directory.
The QC reports are located in the multiqc folder.

### Test data
You may want to test the installation before going with your own data. 

If so, you use use ````--inputcsv tests/fastq/input.csv```` e.g.
````
cd /home/${USER}/work/ssdsnextflowpipeline
conda activate nextflow-dev
sbatch -p computepart -J SSDSnextflowPipeline --export=ALL --mem 5G -t 5-0:0 --mem-per-cpu=1000 --wrap "nextflow run main.nf -c conf/igh.config --inputcsv /home/${USER}/work/ssdsnextflowpipeline/tests/fastq/input.csv --name your_analysis_name --genome mm10 -profile conda"
````
or use ``run_pipeline.sh`` script :
````
bash run_pipeline.sh 
````
**Note : the default value of ``--with_sds_multiqc`` is set to false. If you want to use the SSDS multiQC you need to create a conda environment ; activate the environment, then build the libraries and run the pipeline with ``--with_ssds_multiqc`` and ``--multiqc_dev_conda_env path/to/the/conda/env``** 

