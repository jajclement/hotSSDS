## **SSDS nextflow pipeline version 2.0**
### **Parse, Align and Call Hotspots from single stranded DNA reads from SSDS ChIP-Seq**
**\/!\ Work in progress /!\**


### **Welcome**
**This pipeline is based on  [SSDS pipeline](https://github.com/kevbrick/SSDSnextflowPipeline) and [SSDS call peaks pipeline](https://github.com/kevbrick/callSSDSpeaks) by Kevin Brick, updated and adapted to IGH cluster.**  
See [initial paper](https://genome.cshlp.org/content/22/5/957.long) and [technical paper](https://www.sciencedirect.com/science/article/pii/S0076687917303750?via%3Dihub).  
The pipeline uses [Nextflow]( https://www.nextflow.io/) > 20.04.1  
Briefly, the update from SSDS pipeline version 1.8_NF included **conda profile**, input modification, callpeaks and IDR procedure addition, and global nextflow homogeneisation.  

### **Important :** The pipeline takes **paired-end** reads as input, in fastq(.gz) format.

*Please report any bug occured during the installation to pauline.auffret at igh dot cnrs dot fr*  
  
**The main steps of the pipeline include :**   
* Quality & adapters trimming of paired-end raw fastq files
* Separate mapping of trimmed R1 and R2 reads to the reference genome using bwa & bwa-ra (right align) algorithms
* Bam files filtering
* Parsing of filtered bam files into [!5 main categories](https://genome.cshlp.org/content/22/5/957/F1.large.jpg) : 
* * type 1 reads (T1) : high confidence single stranded DNA
* * type 2 reads (T2) : low confidence single stranded DNA
* * dsDNA : low confidence double stranded DNA
* * dsDNA strict : higher confidence double stranded DNA
* * unclassified : other ambiguous DNA
* Bigwig files generation for visualization through genome browser (for T1 and T2 (optional), options available)
* Peak calling from shuffled T1 bed files 
* Optional IDR (Irreproducible Discovery Rate) analysis to assess replicate consistency) **available for n=2 replicates**
* Peak normalization
* Saturation curve  

In detail, the pipeline is composed of 25 processes :  
* PROCESS 1 : check_design (check input design file)
* PROCESS 2 : makeScreenConfigFile (make configuration file for fastqscreen)
* PROCESS 3 : trimming (use trimmomatic or trim-galore to quality trim, remove adapters and hard trim sequences)
* PROCESS 4 : fastqc (quality control on raw reads using fastqc and fastqscreen)
* PROCESS 5 : bwaAlign (use bwa and custom bwa (bwa right align) to align ssds data)
* PROCESS 6 : filterBam (mark duplicates, remove supplementary alignments, sort and index)
* PROCESS 7 : parseITRs (parse bam files to get the 5 different ssds types : sstype1, sstype2, ds, dsstrict, unclassified)
* PROCESS 8 : makeBigwig (generate normalized bigwig files for T1 and T1+T2 bed files)
* PROCESS 9 : makeBigwigreplicates (optional, generate normalized bigwig files for merged replicates T1 and T1+T2 bed files)
* PROCESS 10 : makeDeeptoolsbigwig (optional, generates bigwig files, coverage and cumulative coverage plots)
* PROCESS 11 : toFRBigwig (optional, generates fwd/rev bigwig files)
* PROCESS 12 : samStats (generates samstats reports)
* PROCESS 13 : makeSSreport (compute stats from ssds parsing)
* PROCESS 14 : ssds_multiqc (make multiqc report for ssds files)
* PROCESS 15 : makeFingerprint (make deeptools fingerprint plots)
* PROCESS 16 : shufBEDs (T1 bed shuffling)
* PROCESS 17 : callPeaks (peak calling on T1 with macs2)
* PROCESS 18 : createPseudoreplicates (optional, creates all pseudoreplicates and pool for idr analysis)
* PROCESS 19 : callPeaksforIDR (optional, call peaks with mac2 on all replicates and pseudo replicates)
* PROCESS 20 : IDRanalysis (optional, perform idr analysis on 4 pairs of replicates or pseudoreplicates
* PROCESS 21 : IDRpostprocess (optional, idr peaks post processing)
* PROCESS 22 : normalizePeaks (center and normalize peaks)
* PROCESS 23 ; mergeFinalPeaks (merge peaks from replicates)
* PROCESS 24 : makeSatcurve (optional, create saturation curve)
* PROCESS 25 : general_multiqc (generates general multiqc report)  

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
sbatch -p computepart -J "install_conda_env" --export=ALL --mem 5G -t 5-0:0 --wrap "bash src/install_pipeline.sh"
```` 
This will create 1 conda environment : **nextflow_dev**.
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
Replicate samples must have the same "group" ID ; the same "antibody" and the same control group.
Control (input) samples must have the 2 last fields ("antibody" and "control") empty. 

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
conda activate nextflow_dev
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


Please use the bash launching script ``run_pipeline.sh`` to execute the pipeline. 
See the available options before use with :
````
bash run_pipeline.sh -h
````
Then run the pipeline with :
````
bash run_pipeline.sh -i inputfile.csv -y "--name your_analysis_name --genome mm10 --profile conda --with_control --nb_replicates 2"
````

**You can use ``-resume`` option (Nextflow native options) to prevent the entire workflow to be rerun in case you need to relaunch an aborted workflow** 




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
A small dataset is present in ````tests/fastq```` directory. 

To use it, you can run :
````
bash run_pipeline.sh -t 1
````
**Without -i option**. Please check -p (pipeline directory) ; -a (conda environment path) and -o (output directory) parameters before run  

Use ````bash run_pipeline.sh -h````

This should take around 10 minutes to run.

### Notes
The default value of ``--with_sds_multiqc`` is set to false. If you want to use the SSDS multiQC you need to create a conda environment ; activate the environment, then build the libraries and run the pipeline with ``--with_ssds_multiqc`` and ``--multiqc_dev_conda_env path/to/the/conda/env``** 

