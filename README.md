# 🐁 **SSDS nextflow pipeline version 2.0** 🐭
**Parse, Align and Call Hotspots from single stranded DNA reads originated from SSDS (Single Stranded Dna Sequencing)**   

## **Context**  
**SSDS method** was originally published by [Khil et al., 2012](https://genome.cshlp.org/content/22/5/957.long).  
The objective is to map **double-strand breaks** (DSBs) along the genome.   
In this method, chromatin is extracted from adult testes and then immunoprecipitated with an antibody against **DMC1 protein**, which is a meiosis-specific recombinase. DMC1 covers the **single-stranded DNA** resulting from the **resection of double-strand breaks (DSBs)**. SSDS uses the ability of single-stranded DNA to form hairpins. The pipeline processes these specific data and identifies **recombination hotspots**.      

See [review](https://pubmed.ncbi.nlm.nih.gov/24136506/) to learn more about meiotic recombination.

     
⚠️ **Work in progress** ⚠️   


*Please give me some feedback if you are using this pipeline, thank you* :blush:    

## **Pipeline overview**  
**This pipeline is based on  [SSDS pipeline](https://github.com/kevbrick/SSDSnextflowPipeline) and [SSDS call peaks pipeline](https://github.com/kevbrick/callSSDSpeaks) by [Kevin Brick, PhD (NIH)](https://www.researchgate.net/profile/Kevin-Brick), updated and adapted to IGH cluster.**  
See [initial paper](https://genome.cshlp.org/content/22/5/957.long) and [technical paper](https://www.sciencedirect.com/science/article/pii/S0076687917303750?via%3Dihub).  
The pipeline uses [Nextflow]( https://www.nextflow.io/) > 20.04.1  
Briefly, the update from SSDS pipeline version 1.8_NF included **conda profile**, input modification, callpeaks, post-processing and IDR procedure addition, and global nextflow homogeneisation.   

  
**The main steps of the pipeline include :**
- Quality check, hard & adapters trimming of paired-end raw fastq files
- Three-step mapping of trimmed R1 and R2 reads to the reference genome :        
	- aligns R1 reads (fully complementary to the genome) with *bwa aln*
	- aligns R2 reads (potentially contain fill-in ITR part at the end of the 5') with *bwa-ra aln* (custom version of bwa that search for the longest mappable suffix in the query)
	- merge with *bwa sampe*
- Bam files filtering for duplicates, unmapped and unpaired reads
- Parsing of filtered bam files into [5 main categories](https://genome.cshlp.org/content/22/5/957/F1.large.jpg) : 
	- type 1 reads (T1) : high confidence single stranded DNA (ITR > 5bp AND microHomology > 0bp AND fill-in > 2bp)
	- type 2 reads (T2) : low confidence single stranded DNA (ITR > 5bp AND microHomology > 0bp AND fill-in < 3bp)
	- dsDNA : low confidence double stranded DNA (ITR < 3bp)
	- dsDNA strict : higher confidence double stranded DNA (ITR < 1bp AND microHomology < 1bp)
	- unclassified : everything that remains unclassified 
- Bigwig files generation for visualization through genome browser (for T1 and T2 (optional), options available)
- Peak calling from shuffled T1 bed files 
- Optional IDR (Irreproducible Discovery Rate) analysis to assess replicate consistency) **available for n=2 replicates**
- Peak normalization
- Saturation curve   
     

**In details, the pipeline is composed of 26 processes :**      
- **Section 1 : INPUT SETTINGS**
	- Process 1  : check_design (check input design file)
	- Process 2  : makeScreenConfigFile (make configuration file for fastqscreen)   
- **Section 2 : TRIMMING**         
	- Process 3  : crop (hard trimming and quality control on raw reads using fastqc and fastqscreen)
	- Process 4  : trimming (use trimmomatic or trim-galore to quality trim and remove adapters from raw sequences)   
- **Section 3 : MAPPING AND PARSING**   
	- Process 5  : bwaAlign (use bwa and custom bwa (bwa right align) to align ssds data)
	- Process 6  : filterBam (mark duplicates, remove supplementary alignments, sort and index)
	- Process 7  : parseITRs (parse bam files to get the different ssds types)
	- Process 8  : makeBigwig (generate normalized bigwig files for t1 and t1+t2 bed files)   
- **Section 4 : BIGWIG**                
	- Process 9  : makeBigwigReplicates (optional, generate normalized bigwig files for merged replicates t1 and t1+t2 bed files)
	- Process 10 : makeDeeptoolsBigWig (optional, generate bigwig files, coverage and cumulative coverage plots)
	- Process 11 : toFRBigWig (optional, generate fwd/rev bigwig files)   
- **Section 5 : PEAK CALLING**          
	- Process 12 : shufBEDs (bed shuffling)
	- Process 13 : callPeaks (peak calling with macs2)   
- **Section 6 : SSDS QC**               
	- Process 14 : samStats (generates samstats reports)
	- Process 15 : makeSSreport (compute stats from ssds parsing)
	- Process 16 : makeFingerPrint (make deeptools fingerprint plots)
	- Process 17 : computeFRIP (Optional, compute frip score for parsed bam file for new genomes)
	- Process 18 : ssds_multiqc (make multiqc report for ssds files)    
- **Section 7 : OPTIONAL IDR ANALYSIS**      
	- Process 19 : createPseudoReplicates (optional, creates all pseudoreplicates and pool for idr analysis)
	- Process 20 : callPeaksForIDR (optional, call peaks with mac2 on all replicates and pseudo replicates)
	- Process 21 : IDRanalysis (optional, perform idr analysis on 4 pairs of replicates or pseudoreplicates)
	- Process 22 : IDRpost- Process (optional, idr peaks postprocessing)     
- **Section 8 : PEAK POST-PROCESSING**     
	- Process 23 : normalizePeaks (center and normalize peaks)
	- Process 24 : mergeFinalPeaks (merge peaks from replicates)
	- Process 25 : makeSatCurve (optional, create saturation curve)  
- **Section 9 : GENERAL QC**            
	- Process 26 : general_multiqc (generate general multiqc report)    

    
## **How to run the pipeline on IGH cluster**
### 0. Requirements
Get [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) if not installed on your system. First, download the installation script (for example in ``/home/${USER}/work/bin`` directory) :
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
First you need to clone the pipeline in your working directory (in the following instructions, ``/home/${USER}/work/`` will refer to your working directory. Please substitute with the according path if different) :
````sh
cd /home/${USER}/work
git clone https://gitlab.igh.cnrs.fr/pauline.auffret/ssdsnextflowpipeline.git
cd ssdsnextflowpipeline
````
Then install pipeline conda environment through sbatch script : 
```` 
sbatch -p computepart -J "install_conda_env" --export=ALL --mem 5G -t 5-0:0 --wrap "bash src/install_pipeline.sh"
```` 
Please use ``bash src/install_pipeline.sh -h`` to see details.   
This will create 1 conda environment : **nextflow_dev** containing **nextflow version 20.04.01**.  
You can use your own conda environment as long as it contains that version of Nextflow.   

### 2. Pipeline configuration 
There are currently **3 configuration files** :
- ````./conf/igh.config```` contains cluster resources requirements & reference genomes info. You don't need to edit this file unless you want to custom the requirements for CPU/memory usage and compute queue (see [IGH cluster documentation](https://kojiki.igh.cnrs.fr/doku.php?id=cluster,)). If a new genome is available on the cluster and does not appear in this file, please contact me or Aubin Thomas, manager of bioinformatics resources at the IGH. If you are running the pipeline on an other computing cluster, you need to specify the relevant configuration file.
- ````./nextflow.config```` contains default pipeline parameters (it's better to **not** edit this file, default parameters will be overwritten by your custom json parameter file, see next point).
- ````./conf/mm10.json```` custom parameters file, you can use your own json file following the same template.    
    
The complete list of parameters is accessible through the command :
````
cd /home/${USER}/work/ssdsnextflowpipeline
conda activate nextflow_dev
nextflow run main.nf --help
````
````
Input data parameters:
    --inputcsv                  FILE    PATH TO INPUT CSV FILE (template and default : /work/demassyie/ssdsnextflowpipeline/tests/fastq/input.csv)
    -params_file                FILE    PATH TO PARAMETERS JSON FILE (template and default : /work/demassyie/ssdsnextflowpipeline/conf/mm10.json)
    --genomebase                DIR     PATH TO REFERENCE GENOMES (default : "/poolzfs/genomes")
    --genome                    STRING  REFERENCE GENOME NAME (must correspond to an existing genome in your config file, default : "mm10")
    --genomedir                 DIR     PATH TO GENOME DIRECTORY (required if your reference genome is not present in your config file)
    --genome_name               STRING  REFERENCE GENOME NAME (e.g ".mm10", required if your reference genome is not present in your config file)
    --genome_fasta              FILE    PATH TO FILE GENOME FASTA FILE WITH PREEXISTING INDEX FILES FOR BWA (required if your reference genome is not present in your config file)
    --fai                       FILE    PATH TO GENOME FAI INDEX FILE (required if your reference genome is not present in your config file)
    --genome2screen             STRING  GENOMES TO SCREEN FOR FASTQC SCREENING (default : ['mm10','hg19','dm3','dm6','hg38','sacCer2','sacCer3'], comma separated list of genomes to screen reads for contamination, names must correspond to existing genomes in your config file)
    --chrsize                   FILE    Chromosome sizes file, default : /work/demassyie/ssdsnextflowpipeline/accessoryFiles/SSDS/mm10/mm10.chrom.sizes (downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes 2021-01-11)
    --hotspots                  DIR     PATH TO HOTSPOTS FILES DIRECTORY (set to "None" if none provided ; default : /work/demassyie/ssdsnextflowpipeline/accessoryFiles/SSDS/hotspots)
    --blacklist                 FILE    PATH TO BLACKLIST BED FILE FOR PEAK CALLING AND IDR (set to "None" if none provided ; default : /work/demassyie/ssdsnextflowpipeline/accessoryFiles/SSDS/blacklist/mm10/blackList.bed)

Output and temporary directory parameters:
    --name                      STRING  ANALYSIS NAME (default : "SSDS_pipeline")
    --outdir                    DIR     PATH TO OUTPUT DIRECTORY (default : name.outdir)
    --publishdir_mode           STRING  MODE FOR EXPORTING PROCESS OUTPUT FILES TO OUTPUT DIRECTORY (default : "copy", must be "symlink", "rellink", "link", "copy", "copyNoFollow","move", see https://www.nextflow.io/docs/latest/process.html)
    --scratch                   DIR     PATH TO TEMPORARY DIRECTORY (default : scratch)

Pipeline dependencies:
    --src                       DIR     PATH TO SOURCE DIRECTORY (default : /work/demassyie/ssdsnextflowpipeline/accessoryFiles/SSDS/scripts ; contains perl scripts)
    --custom_bwa                EXE     PATH TO CUSTOM BWA EXEC (default : /work/demassyie/ssdsnextflowpipeline/accessoryFiles/SSDS/bwa_0.7.12)
    --custom_bwa_ra             EXE     PATH TO CUSTOM BWA_SRA EXEC (default : /work/demassyie/ssdsnextflowpipeline/accessoryFiles/SSDS/bwa_ra_0.7.12)

Trimming parameters:
    --with_trimgalore           BOOL    Use trim-galore instead of Trimmomatic for quality trimming process (default : false)
    --trimgalore_adapters       FILE    trim-galore : PATH TO ADAPTERS FILE (default : none)
    --trimg_quality             INT     trim-galore : minimum quality (default 10)
    --trimg_stringency          INT     trim-galore : trimming stringency (default 6)
    --trim_minlen               INT     trimmomatic : minimum length of reads after trimming (default 25)
    --trim_cropR1               INT     fastx : Cut the R1 read to that specified length (default 50)
    --trim_cropR2               INT     fastx : Cut the R2 read to that specified length (default 50)
    --trim_slidingwin           STRING  trimmomatic : perform a sliding window trimming, cutting once the average quality within the window falls below a threshold (default "4:15")
    --trim_illumina_clip        STRING  trimmomatic : Cut adapter and other illumina-specific sequences from the read (default "2:20:10")
    --trimmomatic_adapters      FILE    PATH TO ADAPTERS FILE FOR TRIMMOMATIC (default /work/demassyie/ssdsnextflowpipeline/TruSeq2-PE.fa, special formatting see http://www.usadellab.org/cms/?page=trimmomatic)

Mapping parameters:
    --with_multimap             BOOL    Keep multimapping reads from bam (default : false)
    --bamPGline                 STRING  bam header (default '@PG        ID:ssDNAPipeline2.0_PAUFFRET')
    --filtering_flag            INT     SAM flag for filtering bam files (default : 2052 ; see https://broadinstitute.github.io/picard/explain-flags.html)
    --picard_min_distance       INT     Picard parameter for marking duplicates (--MINIMUM_DISTANCE) :  width of the window to search for duplicates of a given alignment, default : -1 (twice the first read's read length)
    --picard_optdup_distance    INT     Picard parameter for marking duplicates (--OPTICAL_DUPLICATE_PIXEL_DISTANCE) : The maximum offset between two duplicate clusters in order to consider them optical duplicates. The default is appropriate for unpatterned versions of the Illumina platform (HiSeq2500). For the patterned flowcell models (Novaseq 6000), 2500 is more appropriate, default : 100
    --get_supp                  BOOL    Publish bam files for supplementary aligments, default : false

Bigwig parameter:
    --bigwig_profile            STRING  Bigwig profile using  bedtools (normalization by total library size) : "T1" will produce bigwig for T1 bed files only, one per replicates ; "T12" will also produce bigwig for merged T1+T2, one per replicates ; "T1rep" will also produce T1 bigwig for merged replicates ; "T12rep" will also produce T1+T2 bigwig for merged replicates (default : "T1")
    --kbrick_bigwig             BOOL    Compute bigwig files ; FR bigwig files and coverage plots as in original pipeline by Kevin Brick using deeptools with FPKM normalization (default : false)
    --binsize                   INT     Deeptools binsize parameter (used only if kbrick_bigwig is TRUE ; default : 50)

Peak calling parameters:
    --with_control              BOOL    Use input control files for peak calling analysis (default : false)
    --satcurve                  BOOL    Plot saturation curve (default : false)
    --sctype                    STRING  Saturation curve type (either 'minimal', 'standard' or 'expanded' ; default : 'standard')
    --reps                      INT     Number of iterations for saturation curve (default : 3)
    --bed_trimqual              INT     Mapping quality threshold for bed filtering (default : 30)
    --macs_bw                   INT     Macs2 callpeak bandwidth parameter (default : 1000)
    --macs_slocal               INT     Macs2 callpeak slocal parameter (default : 5000)
    --macs_extsize              INT     Macs2 callpeak extsize parameter (default : 800)
    --macs_qv                   FLOAT   Macs2 callpeak q-value parameter (default : 0.1)
    --macs_pv                   FLOAT   Macs2 callpeak p-value parameter, if not -1, will overrule macs_qv, see macs2 doc (default : -1)
    --no_chrY                   BOOL    Filter out chromosomeY peaks from final peak bed files (default : true)

Optional IDR analysis parameters (ENCODE procedure, see https://github.com/ENCODE-DCC/chip-seq-pipeline2) :
    --with_idr                  BOOL    Perform IDR analysis, only possible if nb_replicates=2 (default : false)
    --nb_replicates             INT     Number of replicates per sample (default : 2)
    --idr_peaktype              STRING  The peak file format for IDR (narrowPeak, regionPeak or broadPeak, default : "regionPeak")
    --idr_threshold             FLOAT   idr p-value threshold (default : 0.1)
    --idr_rank                  INT     p.value or q.value (default : p.value)
    --idr_filtering_pattern     STRING  Regex for filtering bed files (default :"chr[1-9X]+")
    --idr_macs_qv               FLOAT   Macs2 callpeak q-value parameter (default : -1)
    --idr_macs_pv               FLOAT   Macs2 callpeak p-value parameter, if not -1, will overrule macs_qv, see macs2 doc (default : 0.1)

QC parameters:
    --with_ssds_multiqc         BOOL    RUN SSDS MULTIQC (need a functional conda environment, see multiqc-dev_conda-env parameter ; default : false)
    --multiqc_dev_conda_env     DIR     PATH TO MULTIQC-DEV CONDA ENVIRONMENT (used when --with_ssds-multiqc is true ; default : multiqc_dev)
    --multiqc_configfile        FILE    OPTIONAL : PATH TO MULTIQC CUSTOM CONFIG FILE (default : /work/demassyie/ssdsnextflowpipeline/multiqc_config.yaml)

Nextflow Tower parameter:
    -with-tower                 BOOL    Enable job monitoring with Nextflow tower (https://tower.nf/)
    --tower_token               STRING  Nextflow tower key token (see https://tower.nf/ to create your account)

````
     
One **important thing to note**, in Nextflow command lines, the **native options** are preceded with one **single hyphen** (e.g. ``-profile``), while **parameters specific to SSDS pipeline** are preceded with **2 hyphens** (e.g. ``--genome 'mm10'``).   
   

### 3. Input data

1. **Raw reads**     
The pipeline will only process **paired-end data** in fastq(.gz) format.
The input data must be described in an input csv file with the 6 following fields (see template in ``/home/${USER}/work/ssdsnextflowpipeline//tests/fastq/input.csv``)
````
group,replicate,fastq_1,fastq_2,antibody,control
DMC1-chip,1,/work/${USER}/data/SRR1035576_R1.fastq.gz,/work/${USER}/data/SRR1035576_R2.fastq.gz,antiDMC1,Input
DMC1-chip,2,/work/${USER}/data/SRR1035577_R1.fastq.gz,/work/${USER}/data/SRR1035577_R2.fastq.gz,antiDMC1,Input
Input,1,/work/${USER}/data/SRR1035578_R1.fastq.gz,/work/${USER}/data/SRR1035578_R2.fastq.gz,,
````
**Replicate samples** must have the **same "group" ID** ; the **same "antibody"** and the **same control group**.    
Control (input) samples must have the 2 last fields ("antibody" and "control") empty. 

2. **Reference genome**    
The reference genome should be in the ``/poolzvs/genomes`` directory on IGH cluster. Currently, the available genomes are mm10, hg19, hg38, sacCer2, sacCer3, dm3, dm6.

You can use ````--genome 'mm10'```` and you won't need to worry about the other genome parameters like fasta path etc.

But if you want to use **another reference**, you will need to set the following parameters : 
- absolute path to genome ````--genomedir /path/to/genome````
- absolute path to fasta file. **Indexes for BWA SHOULD EXIST in the same directory** ````--genome_fasta /path/to/genome.fa````
- the name of the genome ````--genome_name mm11````
- absolute path to the fai index ````--fai /path/to/genome.fa.fai````
- absolute path the chromosome size file ````--chrsize /path/to/genome.chrom.size````  
     
Do not forget to set accordingly the parameters for ``--hotspots`` and ``--blacklist``.    

### 4. Run the pipeline !
   
There are **2 ways** for running the pipeline.
1. **Run through bash script ``bash run_pipeline.sh -i input_file [options]`` (RECOMMENDED)**    
Options :
	-  ``-h`` display help message    
	- ``-i`` Absolute path to input csv file (REQUIRED except in case of run test on test dataset) **This option matches ``--inputcsv`` parameter in nextflow command line**    
	- ``-g`` Absolute path to the genome config file (default : ``/home/${USER}/work/ssdsnextflowpipeline/conf/mm10.json``) **This option matches ``-params_file`` parameter in nextflow command line**   
	- ``-p`` Absolute path to ssds nextflow pipeline base directory (default : ``/home/${USER}/work/ssdsnextflowpipeline``)    
	- ``-b`` Absolute path to base directory where the output directory will be created (default : ``/home/${USER}/work/results``) **This option matches ``--outdir`` parameter in nextflow command line**   
	- ``-n`` Analysis name (default : SSDS_pipeline) **this option matches ``--name`` parameter in nextflow command line**    
	- ``-c`` Absolute path to IGH cluster configuration file (default :`` /home/${USER}/work/ssdsnextflowpipeline/conf/igh.config``) **This option matches ``-c`` parameter in nextflow command line**     
	- ``-a`` Absolute path to conda environment for nextflow (default : ``/home/${USER}/work/bin/miniconda3/envs/nextflow_dev``)     
	- ``-o`` Optional arguments for the pipeline (for example ``"--with_control --no_multimap --trim_cropR1 50 --trim_cropR2 50"`` ;  default : ``"-profile conda "``) **Be cautious with the ``--`` or ``-``, see section 2.      
	- ``-t`` set to 1 if running pipeline on test data located in ``/home/${USER}/work/ssdsnextflowpipeline/tests/fastq`` (default : 0)      
	- ``-f`` set to 1 to force pipeline to run without checking resume/output directory (default : 0)       
INFO : the output directory will be located in the base directory and will be named after the analysis name parameter with the .outdir suffix (default ``/home/${USER}/work/results/SSDS_pipeline.outdir``)     

    
Please provide all files and directories with **absolute paths**.   
The command line should look like :
 ````
bash run_pipeline.sh -i inputfile.csv -g custom_params.json -n your-analysis-name -a my_own_conda_env -b my-base-directory -o "-profile conda -resume"
````
The results will be located in a directory named after the analyis name (-n argument) suffixed with .outdir, in the base directory.

Example : 
````
bash run_pipeline.sh -i /home/demassyie/work/results/ChIP_Dmc1_cKO_Hells.outdir/infos/input_M_WT_D1_HE_KO.csv -g /home/demassyie/work/results/ChIP_Dmc1_cKO_Hells.outdir/infos/mm10.json -n ChIP_Dmc1_cKO_Hells_newtrim -a /home/demassyie/work/bin/miniconda3/envs/nextflow-dev -o "-profile conda -resume" 
````    

2. **Run directly with Nextflow**     
The main parameters that need to be set are :
	-  ``--inputcsv`` : path to the input csv file, see section 3. This option matches the ``-i`` argument in ``run_pipeline.sh`` script.
	- ``-params_file`` : json file containing list of parameters
	- ``-profile conda`` : Run within conda environment (only available option at the moment)
	- ``-resume`` : Prevent the entire workflow to be rerun in case you need to relaunch an aborted workflow
	- ``--name`` : analysis name, e.g. "SSDS_SRA5678_DMC1". This option matches the ``-n`` argument in ``run_pipeline.sh`` script.
	- ``--outdir`` : path to output directory. This option matches ``base-directory/analysis-name.outdir`` (i.e. defined by ``-b`` and ``-n`` arguments in ``run_pipeline.sh`` script).    
And, if not set in json file :    
	- ``--genome`` : the reference genome, see section 3
	- ``--with_control`` (true/false) : use input control files, see section 3
	- ``--no_multimap`` (true/false) : remove multimappers from bam files
	- ``--nb_replicates`` : number of biological replicates you are running with (maximum 2)
	- ``--bigwig_profile`` : indicates which bigwig to generate (T1 ; T12 ; T1rep or T12rep)
	- ``--with_idr`` (true/false) : run IDR analysis (if nb_replicates=2)
	- ``--satcurve`` (true/false) : plot saturation curve     


### 5. Monitor the pipeline with Nextflow Tower
You can use ``-with-tower`` option to monitor your jobs through [nextflow tower web interface](https://tower.nf/).   
You first need to sign in to get your key, then add it to your parameters with the ``--tower-token 'yourkey'`` option.

### 6. Output
The main output folder is specified through the ````-b```` and ````-n```` arguments in ``run_pipeline.sh`` ; i.e. ``base-directory/analysis-name.outdir``    
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

## Test data
You may want to test the installation before going with your own data. 
A small dataset is present in ````tests/fastq```` directory. 

To use it, you can run :
````
bash run_pipeline.sh -t 1
````
**Without -i option**. Please check ``-p`` (pipeline directory) ; ``-a`` (conda environment path) ; ``-n`` (analysis name) and ``-b`` (base directory for output) parameters before running.      

Use ````bash run_pipeline.sh -h```` to see available options.   

This should take around 10 minutes to run.

## Notes and future developpments
The default value of ``--with_sds_multiqc`` is set to false. If you want to use the SSDS multiQC you need to create a conda environment ; activate the environment, then build the libraries and run the pipeline with ``--with_ssds_multiqc`` and ``--multiqc_dev_conda_env path/to/the/conda/env``** 

See WIP (work in progress) directory.   
- Singularity container
- Jupyter notebooks for beautiful QC

   
🐠🌈
