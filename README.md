# 🐭 **hotSSDS nextflow pipeline version 2.0** :fish:🌈
**Parse, Align and Call Hotspots from single stranded DNA reads originated from SSDS (Single Stranded Dna Sequencing)**   

## **Context**  
**SSDS method** was originally published by [Khil et al., 2012](https://genome.cshlp.org/content/22/5/957.long).  
The objective is to map **double-strand breaks** (DSBs) along the genome.   
In this method, chromatin is extracted from adult testes and then immunoprecipitated with an antibody against **DMC1 protein**, which is a meiosis-specific recombinase. DMC1 covers the **single-stranded DNA** resulting from the **resection of double-strand breaks (DSBs)**. SSDS uses the ability of single-stranded DNA to form hairpins. The pipeline processes these specific data and identifies **recombination hotspots**.      

See [this review](https://pubmed.ncbi.nlm.nih.gov/24136506/) by 
Frédéric Baudat, Yukiko Imai and Bernard de Massy to learn more about meiotic recombination.

## **Pipeline overview**  
**This pipeline is based on  [SSDS pipeline](https://github.com/kevbrick/SSDSnextflowPipeline) and [SSDS call peaks pipeline](https://github.com/kevbrick/callSSDSpeaks) by [Kevin Brick, PhD (NIH)](https://www.researchgate.net/profile/Kevin-Brick), updated and adapted to IGH cluster.**  
See [initial paper](https://genome.cshlp.org/content/22/5/957.long) and [technical paper](https://www.sciencedirect.com/science/article/pii/S0076687917303750?via%3Dihub).  
The pipeline uses [Nextflow]( https://www.nextflow.io/) > 20.04.1  
Briefly, the update from SSDS pipeline version 1.8_NF included **conda/mamba/singularity/docker execution profiles**, input modification, callpeaks, post-processing and IDR procedure addition, and global nextflow homogeneisation.   

  
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
- Optional IDR (Irreproducible Discovery Rate) analysis to assess replicate consistency) **only available for n=2 replicates**
- Peak centering and normalization
- Saturation curve   
- Quality control reports generation
     

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

    
## **Requirements**

The hotSSDS nextflow pipeline requires [Nextflow DSL1 version 20.10.0](https://github.com/nextflow-io/nextflow/),  and at least one of the following : 
- [conda](https://docs.conda.io/en/latest/)
- [mamba](https://github.com/mamba-org/mamba)
- [Docker](https://www.docker.com/)
- [Singularity version 3.4.1](https://sylabs.io/)
   
The preferable execution profile is **Singularity** to ensure full portability.

Nextflow can be easily installed using [conda package manager](https://anaconda.org/bioconda/nextflow).

Downloading the pipeline will preferably require [git](https://git-scm.com/), otherwise the pipeline can be downloaded using command-line programs for retrieving files from Internet such as [wget](https://www.gnu.org/software/wget/). Nextflow can easily be installed using conda package manager [https://anaconda.org/bioconda/nextflow].
The hotSSDS pipeline will rely on the pre-existence of genome references on your system. If you need to download them, then [bwa](https://bio-bwa.sourceforge.net/bwa.shtml) and [samtools](http://www.htslib.org/) will be required as well.

### **Download the pipeline**
Using git (recommended)
````sh
git clone https://github.com/jajclement/hotSSDS.git
cd hotSSDS
````
Or download zip file using wget, then unzip file
````sh
wget https://github.com/jajclement/hotSSDS/archive/refs/heads/master.zip
unzip hotSSDS-master.zip
mv hotSSDS-master hotSSDS
cd hotSSDS
````

### **Download Singularity images**
Singularity images are used to encapsulate all required softwares and dependencies for the different steps of the pipeline. They make the pipeline portable on different systems. As they can be voluminous, they are not included in the pipeline git repository.  
  
Prior to run the pipeline using Singularity execution profile, it is necessary to download the images from [Zenodo ‘hotSSDS Pipeline Singularity Images’ open repository](https://zenodo.org/record/7783473).   
  
To do this, two options are available depending on whether the computing environment on which the pipeline will be executed has access to the Internet (see point a) or not (see point b).  

#### a. Run the pipeline using the option ```--get_sif```
The option ```--get_sif``` allows to launch a « dry run » that will check the existence of Singularity images in the pipeline directory. If not present, the pipeline will download them. Once download is completed, the pipeline stops. It can then be run without the option ```--get_sif``` to perfom hotSSDS analyses.
````sh
nextflow run main.nf -c conf/cluster.config \
	-params-file conf/test.json \
	–profile test,<singularity|mamba|conda|docker> \
	--get_sif >& get_sif_main_log.txt 2>&1
````
#### b. Download all singularity images independantly
Download all the 10 *.sif* files from zenodo open repository at https://zenodo.org/record/7783473 and place them in ``hotSSDS/containers`` folder so that the final repository structure is :
````sh
containers/
├── bam-box-1.0
│   ├── bam-box_1.0.sif
│   ├── Dockerfile
│   └── environment.yml
├── bigwig-box-1.0
│   ├── bigwig-box-1.0.sif
│   ├── Dockerfile
│   └── environment.yml
├── frip-box-1.0
│   ├── Dockerfile
│   ├── environment.yml
│   └── frip-box_1.0.sif
├── idr-box-1.0
│   ├── Dockerfile
│   ├── environment.yml
│   └── idr-box_1.0.sif
├── multiqc-box-1.0
│   ├── Dockerfile
│   ├── environment.yml
│   └── multiqc-box_1.0.sif
├── peak-calling-box-1.0
│   ├── Dockerfile
│   ├── environment.yml
│   └── peak-calling-box_1.0.sif
├── plot-box-1.0
│   ├── Dockerfile
│   ├── environment.yml
│   └── plot-box_1.0.sif
├── python-3.8
│   ├── environment.yml
│   └── python-3.8.sif
├── ssds-qc-box-1.0
│   ├── Dockerfile
│   ├── environment.yml
│   └── ssds-qc-box_1.0.sif
└── trimming-box-1.0
    ├── Dockerfile
    ├── environment.yml
    └── trimming-box_1.0.sif
````
### Configure computing parameter files
Edit ``hotSSDS/conf/cluster.config`` file to adjust the parameters to a computing cluster.  
The following sections are expected to be overwritten : 

- **DEFAULT CLUSTER CONFIGURATION** section :  
	- _Cluster description_ subsection (optional) : fill-in cluster description, mail and url ;
    - _Executor (cluster scheduler) parameters_ subsection (mandatory) : in ‘name’ field, specify ‘slurm’ or ‘pbspro’ for instance, depending on your system scheduling program ;
    - _Computing queues parameters_ (mandatory) : in ‘queue’ field, specify the computing queues name(s) on which jobs will be submitted. ;
    - _Max resources_ (mandatory) : adapt max_memory, max_cpus and max_time fields if needed.  


- **PROFILES SPECIFIC PARAMETERS** section :
  - Specify custom launching commands for Conda, Mamba, Singularity and Docker  


- **GENOMES LOCATION** section :  
Write the absolute paths to reference genome(s) :
	- _genome_fasta_ : absolute path to genome fasta file. If you need to, download the fasta file of your reference genome (you can use the [golden path](https://hgdownload.soe.ucsc.edu/downloads.html), the best way to get all the required infos for your genome if it is published on UCSC). There are several ways to do it. One is to directly upload using the wget command (other ways to do it will be specified there later). You often have the instruction to download and even the command line on the webpage where you access the genome.; 
    - _genomedir_ : absolute path to directory containing genome fasta file AND bwa index files (to create using bwa-index) bwa index ref.fa ;
    - _genome_name_	: genome name ;
    - _fai_ : absolute path to genome fai index (to create with [faidx tool](http://www.htslib.org/doc/samtools-faidx.html) :  
  
  	  ``samtools faidx <ref.fasta> -o <ref.fai>``

If needed, edit ``hotSSDS/conf/resources.config``  to adjust specific process resources. To do so, edit cpus/memory/time in **PROCESSES SPECIFIC RESSOURCES REQUIREMENTS** section. You can also add specific computing queues to some categories.

It is important to note that many institutes have one such configuration file referenced in [nf-core/configs repository](https://github.com/nf-core/configs) that you can download and adapt to hotSSDS pipeline.

### **Prepare input file**

The pipeline will only process **paired-end data** in fastq(.gz) format.
The input data must be described in an input csv file with the 6 following fields :
````
group,replicate,fastq_1,fastq_2,antibody,control
DMC1-chip-WT,1,/path/to/data/SRR1035576_R1.fastq.gz,/path/to/data/SRR1035576_R2.fastq.gz,antiDMC1,Input-WT
DMC1-chip-WT,2,/path/to/data/SRR1035577_R1.fastq.gz,/path/to/data/SRR1035577_R2.fastq.gz,antiDMC1,Input-WT
DMC1-chip-KO,1,/path/to/data/SRR1035578_R1.fastq.gz,/path/to/data/SRR1035578_R2.fastq.gz,antiDMC1,Input-KO
DMC1-chip-KO,2,/path/to/data/SRR1035579_R1.fastq.gz,/path/to/data/SRR1035579_R2.fastq.gz,antiDMC1,Input-KO
Input-WT,1,/path/to/data/SRR1035580_R1.fastq.gz,/path/to/data/SRR1035580_R2.fastq.gz,,
Input-KO,1,/path/to/data/SRR1035581_R1.fastq.gz,/path/to/data/SRR1035581_R2.fastq.gz,,
````
**Replicate samples** must have the **same "group" ID** ; the **same "antibody"** and the **same control group**.   
If the samples do not have an associated input control sample, leave the last fields empty.  
Control (input) samples must have the 2 last fields ("antibody" and "control") empty.  
**There must be no empty line at the end of the csv file**

### **Edit/create parameter file**
The complete list of parameters is accessible through the command :
````
nextflow run main.nf --help
````
````
N E X T F L O W  ~  version 20.10.0
Launching `main.nf` [reverent_perlman] - revision: 5244f2a6f5
=============================================================================
  hotSSDS pipeline version 2.0 : Align, parse and call hotspots from SSDNA
=============================================================================
    Usage:

    nextflow run main.nf -c conf/cluster.config --params_file conf/mm10.json --inputcsv tests/fastq/input.csv  --name "runtest" --trim_cropR1 36 --trim_cropR2 40 --with_trimgalore -profile singularity -resume

    Runs with Nextflow DSL1 v20.10.0
=============================================================================
Input data parameters:
    --inputcsv                  FILE    PATH TO INPUT CSV FILE (template and default : hotSSDS/tests/fastq/input.csv)
    -params_file                FILE    PATH TO PARAMETERS JSON FILE (template and default : hotSSDS/conf/mm10.json)
    --genomebase                DIR     PATH TO REFERENCE GENOMES
    --genome                    STRING  REFERENCE GENOME NAME (must correspond to an existing genome in your config file, default : "mm10")
    --genomedir                 DIR     PATH TO GENOME DIRECTORY (required if your reference genome is not present in your config file)
    --genome_name               STRING  REFERENCE GENOME NAME (e.g ".mm10", required if your reference genome is not present in your config file)
    --genome_fasta              FILE    PATH TO GENOME FASTA FILE WITH PREEXISTING INDEX FILES FOR BWA (required if your reference genome is not present in your config file)
    --fai                       FILE    PATH TO GENOME FAI INDEX FILE (required if your reference genome is not present in your config file)
    --genome2screen             STRING  GENOMES TO SCREEN FOR FASTQC SCREENING (default : ['mm10','hg19','dm3','dm6','hg38','sacCer2','sacCer3'], comma separated list of genomes to screen reads for contamination, names must correspond to existing genomes in your config file)
    --chrsize                   FILE    Chromosome sizes file, default : ssdsnextflowpipeline/data/mm10/mm10.chrom.sizes (downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes 2021-01-11)
    --hotspots                  DIR     PATH TO HOTSPOTS FILES DIRECTORY (set to "None" if none provided ; default :  hotSSDS/data/hotspots/mm10/hotspots)
    --blacklist                 FILE    PATH TO BLACKLIST BED FILE FOR PEAK CALLING AND IDR (set to "None" if none provided ; default : hotSSDS/data/blacklist/mm10/blackList.bed)

Output and temporary directory parameters:
    --name                      STRING  ANALYSIS NAME (default : "hotSSDSPipeline")
    --outdir                    DIR     PATH TO OUTPUT DIRECTORY (default : hotSSDS/{params.name}.outdir/02_results")
    --publishdir_mode           STRING  MODE FOR EXPORTING PROCESS OUTPUT FILES TO OUTPUT DIRECTORY (default : "copy", must be "symlink", "rellink", "link", "copy", "copyNoFollow","move", see https://www.nextflow.io/docs/latest/process.html)

Pipeline dependencies:
    --src                       DIR     PATH TO SOURCE DIRECTORY (default : hotSSDS/bin ; contains perl scripts)
    --custom_bwa                EXE     PATH TO CUSTOM BWA EXEC (default : hotSSDS/bin/bwa_0.7.12)
    --custom_bwa_ra             EXE     PATH TO CUSTOM BWA_SRA EXEC (default : hotSSDS/bin/bwa_ra_0.7.12)

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
    --trimmomatic_adapters      FILE    PATH TO ADAPTERS FILE FOR TRIMMOMATIC (default hotSSDS/data/TruSeq2-PE.fa, special formatting see http://www.usadellab.org/cms/?page=trimmomatic)

Mapping parameters:
    --with_multimap             BOOL    Keep multimapping reads from bam (default : false)
    --bamPGline                 STRING  bam header (default '@PG\tID:ssDNAPipeline2.0_PAUFFRET')
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
    --idr_setup                 STRING  Threshold profile for idr. This will define the thresholds for true replicates, pool replicates, self replicates r1 and self replicates r2. Profile "auto" is based on ENCODE guidelines and profile "custom" allows to set custom thresholds (see parameters --idr_threshold_r1 --idr_threshold_r2 --idr_threshold_truerep and --idr_threshold_poolrep ; default : auto)
    --idr_threshold_r1          FLOAT   idr threshold for self replicates r1 (used if --idr_setup is "custom" only ; default : 0.05)
    --idr_threshold_r2          FLOAT   idr threshold for self replicates r2 (used if --idr_setup is "custom" only ; default : 0.05)
    --idr_threshold_truerep     FLOAT   idr threshold for true replicates (used if --idr_setup is "custom" only ; default : 0.05)
    --idr_threshold_poolrep     FLOAT   idr threshold for pooled replicates (used if --idr_setup is "custom" only ; default : 0.01)
    --idr_rank                  INT     p.value or q.value (default : p.value)
    --idr_filtering_pattern     STRING  Regex for filtering bed files (default :"chr[1-9X]+" for mouse ; set ".*" to keep everything)
    --idr_macs_qv               FLOAT   Macs2 callpeak q-value parameter (default : -1)
    --idr_macs_pv               FLOAT   Macs2 callpeak p-value parameter, if not -1, will overrule macs_qv, see macs2 doc (default : 0.1)

QC parameters:
    --with_ssds_multiqc         BOOL    RUN SSDS MULTIQC (default : true)
    --multiqc_configfile        FILE    OPTIONAL : PATH TO MULTIQC CUSTOM CONFIG FILE (default : hotSSDS/conf/multiqc_config.yaml)

Nextflow Tower parameter:
    -with-tower                 BOOL    Enable job monitoring with Nextflow tower (https://tower.nf/)

Singularity images parameters:
    --get_sif                   BOOL    [REQUIRE INTERNET ACCESS] Check and download singularity images if necessary (if true, pipeline will stops after download. Once downloading has been done, relaunch pipeline with false ; default: false)
    --url_sif                   URL     URL TO PUBLIC SINGULARITY IMAGES REPOSITORY (default : https://zenodo.org/record/7783473/files)

````

Pipeline parameters can be set in two different ways :
- via parameter file, 
- directly in the nextflow command-line (parameters passed though command line will overwrite those in parameter file).  


For _Mus musculus_ based analysis, parameter file ``hotSSDS/conf/mm10.json`` contains default parameters that can be overwritten

One **important thing to note**, in Nextflow command lines, the **native options** are preceded with one **single dash** (e.g. ``-profile``), while **parameters specific to SSDS pipeline** are preceded with **2 dashes** (e.g. ``--genome 'mm10'``).   
   
### **Launch the pipeline**

Once you have set computing config file and arameter file, you can launch the pipeline using the following command-line :  
````sh
nextflow run main.nf -c conf/cluster.config \
	-params-file conf/mm10.json \
	--inputcsv /path/to/input.csv \
	-profile <singularity|mamba|conda|docker> \
	--name "My_workflow_name" >& main_log.txt 2>&1
````
It is highly recommended to launch the command in a batch job on the computing cluster, as its execution will take time and computing resources. It is also recommanded to redirect the output of this main nextflow command-line to an identified log file, which will be usefull to monitor the pipeline execution.

The main parameters that need to be set are :
- ``--inputcsv`` : path to the input csv file
- ``-params_file`` : json file containing list of parameters
- ``-profile`` : <conda|mamba|docker|singularity>
- ``-resume`` : Prevent the entire workflow to be rerun in case you need to relaunch an aborted workflow.
- ``--name`` : analysis name, e.g. "SSDS_SRA5678_DMC1".
- ``--outdir`` : path to output directory.    
And, if not set in json file :    
- ``--genome`` : the reference genome
- ``--with_control`` (true/false) : use input control files
- ``--no_multimap`` (true/false) : remove multimappers from bam files
- ``--nb_replicates`` : number of biological replicates you are running with (maximum 2)
- ``--bigwig_profile`` : indicates which bigwig to generate (T1 ; T12 ; T1rep or T12rep)
- ``--with_idr`` (true/false) : run IDR analysis (if nb_replicates=2)
- ``--satcurve`` (true/false) : plot saturation curve     
- ``--with_sds_multiqc`` (true/false) : generate ssds qc plots

### **Run a short test** 
A small dataset can be used to test if the pipeline is correctly running on your system.  
To do so, run :
````sh
nextflow run main.nf -c conf/cluster.config \
	-params-file conf/test.json \
	–profile test,<singularity|mamba|conda|docker> >& test_main_log.txt 2>&1
````
This test run should approximately take 5 minutes to complete.  

On completion, the end of main log test_main_log.txt should look like :
````
executor >  pbspro (25)
[aa/ff63f4] process > check_design (input.csv)			[100%] 1 of 1 ✔
[1f/c81a5f] process > makeScreenConfigFile (TEST_SSDS)		[100%] 1 of 1 ✔
[f1/341fad] process > crop (TEST_IP_R1_T1)			[100%] 1 of 1 ✔
[fe/67d328] process > trimming (TEST_IP_R1_T1)			[100%] 1 of 1 ✔
[cb/a77283] process > bwaAlign (TEST_IP_R1_T1)			[100%] 1 of 1 ✔
[38/e2e8c6] process > filterBam (TEST_IP_R1_T1)			[100%] 1 of 1 ✔
[af/05397c] process > parseITRs (TEST_IP_R1_T1)			[100%] 1 of 1 ✔
[d1/b74fbc] process > makeBigwig (TEST_IP_R1_T1)		[100%] 1 of 1 ✔
[1f/256f21] process > shufBEDs (TEST_IP_R1)			[100%] 1 of 1 ✔
[3c/f2b936] process > callPeaks (TEST_IP_R1)			[100%] 5 of 5 ✔
[43/655990] process > samStats (TEST_IP_R1_T1)			[100%] 5 of 5 ✔
[18/200309] process > makeSSreport (TEST_IP_R1_T1)		[100%] 1 of 1 ✔
[8f/19ef51] process > makeFingerPrint (TEST_SSDS)		[100%] 1 of 1 ✔
[c8/a77ede] process > ssds_multiqc (TEST_IP_R1_T1)		[100%] 1 of 1 ✔
[85/6daa72] process > normalizePeaks (TEST_IP_R1)		[100%] 1 of 1 ✔
[6a/122837] process > makeSatCurve (TEST_SSDS)			[100%] 1 of 1 ✔
[a7/a00c9c] process > general_multiqc (TEST_SSDS)		[100%] 1 of 1 ✔
Completed at: 10-Mar-2023 13:37:48
Duration    : 3m 57s
CPU hours   : 0.3
Succeeded   : 25
````

### **Monitor the pipeline**
- Using **Nextflow Tower** :
You can use ``-with-tower`` option to monitor your jobs through [nextflow tower web interface](https://tower.nf/).   
You first need to sign in to generate an access token then export  **tower_access_token** environment variable in your computing environment.
- Checking the main log file (main_log.txt in the example command line) using :
````sh
tail –d main_log.txt
executor >  pbspro (1)
[21/3786f7] process > check_design (Nore_input_fi... [100%] 1 of 1 ✔
[57/bc36de] process > makeScreenConfigFile (SSDS_... [100%] 1 of 1 ✔
[d3/c42056] process > crop (WT_R2_T1)                [100%] 4 of 4 ✔
[7d/c06fe8] process > trimming (WT_R2_T1)            [100%] 4 of 4 ✔
[ca/74adfd] process > bwaAlign (WT_R2_T1)            [100%] 4 of 4 ✔
[f0/dec2f6] process > filterBam (WT_R2_T1)           [100%] 4 of 4 ✔
[aa/39f5f4] process > parseITRs (WT_R2_T1)           [100%] 4 of 4 ✔
[79/c5fd62] process > makeBigwig (WT_R2_T1)          [100%] 4 of 4 ✔
[22/a02566] process > makeDeeptoolsBigWig (WT_R2_T1) [100%] 20 of 20 ✔
[1b/218fa3] process > toFRBigWig (WT_R2_T1)          [100%] 20 of 20 ✔
[24/35b6ac] process > shufBEDs (WT_R1)               [100%] 4 of 4 ✔
[78/1c307b] process > callPeaks (MUT_R1)             [100%] 9 of 9 ✔
[18/b8edd6] process > samStats (WT_R2_T1)            [100%] 20 of 20 ✔
[e4/4461d1] process > makeSSreport (WT_R2_T1)        [100%] 4 of 4  ✔
[4b/e769b4] process > makeFingerPrint (SSDS_pipel... [100%] 1 of 1  ✔
[58/02be5a] process > ssds_multiqc (WT_R2_T1)        [100%] 4 of 4  ✔
[5b/be1b53] process > createPseudoReplicates (MUT)   [ 50%] 1 of 2
[37/0acb1c] process > callPeaksForIDR (WT)           [100%] 1 of 1 ✔
[-        ] process > IDRanalysis                    -
[-        ] process > IDRpostprocess                 -
[-        ] process > normalizePeaks_idr             -
[-        ] process > makeSatCurve                   -
[-        ] process > general_multiqc                
````
### Output
The main output folder is specified using ``--outdir`` parameter.   
This folder will contain the following directories :
- **00_reports** : contains Nextflow execution reports and diagram
- **01_logs** : contains log files for all processes 
- **02_results** : contains results files, including :
  * **trimming** with trimmed fastq files
  * **bwa** with mapped reads, filtered mapped reads and parsed fragments in bam (bai) and bed format
  * **bigwig**  with bigwig and bedgraph files
  * **peaks** with raw, filtered, centered peaks in bed format and with saturation curve and IDR files is options are selected
  * **qc** with quality controls files and reports for sequencing, mapping, peak calling and parsing steps
  * **idr** with peak set files for IDR statistical testing
  * **work** is the working directory for Nextflow   

Tree overview of the output folder composition [DEPRECATED] :

````
.
├── bigwig                                              : contains bigwig files according to the parameter set
│   └── T1
│       └── log
├── qc                                                  : contains quality control files, pictures and reports
│   ├── multiqc                                         : contains a summary of QC stats for all processes
│   │   ├── *.multiQC.quality-control.report_plots
│   │   │   ├── svg
│   │   │   ├── png
│   │   │   └── pdf
│   │   └── *.multiQC.quality-control.report_data
│   ├── samstats                                        : contains mapping statistics tabs
│   ├── ssds                                            : contains tabs and plots about SSDS parsing statistics
│   ├── flagstat                                        : contains mapping statistics files
│   ├── fingerprint                                     : contains fingerprint plots
│   ├── trim_fastqc                                     : contains fastqc reports for trimmed reads
│   ├── raw_fastqc                                      : contains fastqc reports for raw reads
│   ├── fastqscreen                                     : contains plots for fastqscreen screening
│   └── design                                          : contains info about the run
│       └── pipeline_info
├── peaks                                               : contains bed files for peaks
│   └── with[out]-input
│       ├── normalized                                  : contains normalized and recentered peaks, 
│       │   └── [no-]idr
│       │       ├── tab
│       │       └── log
│       ├── finalpeaks                                  : contains a copy of final peaks (generally after IDR or merge)
│       ├── saturation_curve                            : contains saturation curve files and plots
│       │   └── standard
│       │       └── peaks
│       ├── macs2                                       : contains raw peaks called by macs2
│       │   └── pv*_qv*_bw*_sloc*_extsize*
│       │       ├── log
│       │       ├── xls
│       │       ├── narrowPeak
│       │       └── bed
│       └── bed_shuffle                                 : contains shuffled bed files before peak calling
│           └── trim_q*
├── bwa                                                 : contains raw, filtered and parsed reads in bam and bed format
│   ├── filterbam
│   │   └── flag_*
│   │       ├── parse_itr
│   │       │   ├── unclassified
│   │       │   │   ├── bed
│   │       │   │   └── bam
│   │       │   ├── type2
│   │       │   │   ├── bed
│   │       │   │   └── bam
│   │       │   ├── type1
│   │       │   │   ├── bed
│   │       │   │   └── bam
│   │       │   ├── norm_factors
│   │       │   ├── log
│   │       │   ├── flagstat
│   │       │   └── dsDNA
│   │       │       └── bed
│   │       ├── log
│   │       ├── bed
│   │       └── bam
│   └── bam
│       └── log
├── trimming                                            : contains trimmed fastq files
│   └── trim_fastq
├── idr                                                 : contains files for IDR process
   └── with[out]-input
       └── narrowPeak_macs2pv*_macs2qv*1_idr_setup-*
           ├── bfilt                                   : contains blacklist filtered bed files
           ├── log
           ├── macs2                                   : contains macs2 narrowpeaks files
           │   └── log
           ├── peaks                                   : contains bed files
           ├── plot                                    : contains IDR plots
           ├── pseudo_replicates                       : contains bed for pseudo replicates files
           ├── qc
           │   └── log
           └── unthresholded-peaks                     : contains unthresholded bed files

````


I recommend to look at ``qc/multiqc/*.multiQC.quality-control.report.html`` file first to have a look at **sequencing, mapping, parsing quality.**   
Then you use the **bigwig files** in your favorite brower or online [IGV](https://igv.org/app/).    
You can also look at the **peaks** in ``peaks/with[out]-input/finalpeaks`` .   
Then, you can run [ssdspostprocess pipeline](https://github.com/jajclement/ssdspostprocess) to go deeper in the peaks analysis.  


### **Debug**
In an ideal world the pipeline would never crash but let's face it, it will happen. 
Here are some clues to help you debug. 
Fisrt, the main log file will give you many precious clues :
- The second line of this file looks like ``Launching `main.nf` [fervent_ptolemy] - revision: e2a189c33c`` indicating that the current run has been internally named **fervent_ptolemy**
- Just after you will find all the parameters used in the current run
- Then you will have the detailled pipeline status, for each process, for example :
````
[88/c511ac] process > check_design (input.csv)       [100%] 1 of 1 ✔
[d7/59a78c] process > makeScreenConfigFile (SSDS_... [100%] 1 of 1 ✔
[75/9a5a8a] process > crop (TEST_IP_R1_T1)           [100%] 1 of 1 ✔
[a7/2b8da0] process > trimming (TEST_IP_R1_T1)       [100%] 1 of 1 ✔
[6c/b564b0] process > bwaAlign (TEST_IP_R1_T1)       [100%] 1 of 1 ✔
[4d/cadedc] process > filterBam (TEST_IP_R1_T1)      [100%] 1 of 1 ✔
[bf/9b1f6a] process > parseITRs (TEST_IP_R1_T1)      [100%] 1 of 1 ✔
[55/89204f] process > makeBigwig (TEST_IP_R1_T1)     [100%] 1 of 1 ✔
[49/a54378] process > shufBEDs (TEST_IP_R1)          [100%] 1 of 1 ✔
[5c/ab370e] process > callPeaks (TEST_IP_R1)         [100%] 5 of 5 ✔
[bc/5e7a98] process > samStats (TEST_IP_R1_T1)       [100%] 5 of 5 ✔
[da/02474b] process > makeSSreport (TEST_IP_R1_T1)   [100%] 1 of 1 ✔
[d0/be7efa] process > makeFingerPrint (SSDS_pipel... [100%] 1 of 1 ✔
[c5/c29f40] process > normalizePeaks (TEST_IP_R1)    [100%] 1 of 1 ✔
[77/ceb36f] process > makeSatCurve (SSDS_pipeline... [100%] 1 of 1 ✔
[72/d60ad9] process > general_multiqc (SSDS_pipel... [100%] 1 of 1 ✔
````
For example, for process **trimming**, the associated key is **a7/2b8da0** meaning that in the nextflow **work** directory, there will be a folder named **a7**, which will contain a folder beginning with **2b8da0** : this is where the detailled log files for process trimming will be. In every log folder like this, you can ckeck the following files (remember to write ``ls -A`` to display files whose names begin with ``.``) :    
    * ``.command.out`` : output of the process    
    * ``.command.err`` : error returned by the process, if any    
    * ``.command.log`` : log file of the process    
    * ``.command.sh`` : bashscript executed for the process    
    * ``.command.run`` : nextflow script executed for the process    
    * ``.command.trace`` : resources used by the process    
    * ``.exitcode`` : exit code of the process (if process succedeed : must be 0)    
    * Optional process specific log (``*.log`` file)
- There is also a ``.nexflow.log`` file created each time the pipeline is run (the latest is named ``.nexflow.log``, the second latest is renamed ``.nexflow.log2`` and so on). This log file can give insights of nextflow - slurm communication during the pipeline, such as jobs ID, run time and so on. This file is located in the folder from where the pipeline is launched.   
- If the pipeline crashed because of time or memory limits, you can edit ``conf/igh.config`` file (usually, exit code is 143)    
- The main causes of pipeline crashed are :    
    * Wrong parameters
    * iles not found
    * Conda drama (if one tool is not found, that means there is an error with conda. Check that pipeline runs using ``-profile conda``
    * Time/memory limits
    * Exotic parameters unforeseen behavior (if running on a new genome for instance)  

             
**Do not hesitate to contact me or open an issue if you can't resolve one bug.**

### **Notes and future developpments**
See [TODO.md](https://github.com/jajclement/hotSSDS/blob/master/TODO.md) file.   


   
:santa:
