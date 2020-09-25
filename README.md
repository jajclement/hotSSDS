## **SSDS nextflow pipeline version 1.8_NF_pa**
**\/!\ Work in progress /!\\**
### **Welcome**
**This is [SSDS pipeline by Kevin Brick](https://github.com/kevbrick/SSDSnextflowPipeline) updated and adapted to IGH cluster.**
See [initial paper](https://genome.cshlp.org/content/22/5/957.long) and [technical paper](https://www.sciencedirect.com/science/article/pii/S0076687917303750?via%3Dihub).
The pipeline uses [Nextflow]( https://www.nextflow.io/) > 20.04.1
Briefly, the update from SSDS pipeline version 1.8_NF included **conda profile**, input and trimming process modification, and global nextflow homogeneisation.

*Please report any bug occured during the installation to pauline.auffret at igh dot cnrs dot fr*


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

### 3. Input data
The pipeline will only process **paired-end data**.
Currently 3 input data formats are supported : 
- fastq(.gz) ````--fqdir input_data/raw_data/*{R1,R2}.fastq.gz````
- bam ````--bamdir input_data/raw_data/*.bam````
- SRA identifiers ````--sra_ids=['ERR908507', 'ERR908506']```` If you use this option, please consider [getting NCBI API key](https://www.ncbi.nlm.nih.gov/myncbi/) and add it to the parameters with ``--ncbi_api_key 'yourapikey'`` **currently not working due to fromSRA bug**

The reference genome should be in the ``/poolzvs/genomes`` directory on IGH cluster. Currently, the available genomes are mm10, hg19, hg38, sacCer2, sacCer3, dm3, dm6.

You can use ````-- genome 'mm10'```` and you won't need to worry about the other genome parameters like fasta path etc.

But if you want to use another reference, you will need to set : 
- path to genome ````--genomedir /path/to/genome````
- path to fasta file. Indexes for BWA SHOULD EXIST in the same directory ````--genome_fasta /path/to/genome.fa````
- the name of the genome ````--genome_name mm11````
- the path to the fai index ````--fai /path/to/genome.fai.fai````

### 4. Run the pipeline !
To see all the parameters, please run :
````
cd /home/${USER}/work/ssdsnextflowpipeline
conda activate nextflow-dev
nextflow run main.nf --help
````
You can either run the pipeline directly through the command line :
````
cd /home/${USER}/work/ssdsnextflowpipeline
conda activate nextflow-dev
sbatch -p computepart -J SSDSnextflowPiepline --export=ALL --mem 5G -t 5-0:0 --mem-per-cpu=1000 --wrap "nextflow run main.nf -c conf/igh.config --fqdir path/to/your/data/*{R1,R2}.fastq.gz --name your_analysis_name --genome mm10 --profile conda
````
or use ``run_pipeline.sh`` script : in that case make sure to edit with your parameters, then
````
bash run_pipeline.sh
````
**Nextflow Tower**
You can use ``-with-tower`` option to monitor your jobs through [nextflow tower web interface](https://tower.nf/). 
You first need to sign in to get your key, then add it to your parameters with the ``--tower-token 'yourkey'`` option.

### Test data
You may want to test the installation before going with your own data. 

If so, just use ````--fqdir tests/fastq/*{R1,R2}.fastq"```` e.g.
````
cd /home/${USER}/work/ssdsnextflowpipeline
conda activate nextflow-dev
sbatch -p computepart -J SSDSnextflowPiepline --export=ALL --mem 5G -t 5-0:0 --mem-per-cpu=1000 --wrap "nextflow run main.nf -c conf/igh.config --fqdir tests/fastq/*{R1,R2}.fastq --name your_analysis_name --genome mm10 -profile conda"
````
or use ``run_pipeline.sh`` script (make sur to edit it with the relevant paramters before use) :
````
bash run_pipeline.sh
````
**Note : the default value of ``--with-ssds_multiqc`` is set to false. If you want to use the SSDS multiQC you need to create a conda environment ; activate the environment, then build the libraries and run the pipeline with ``--with-ssds_multiqc`` and ``--multiqc_dev_conda_env path/to/the/conda/env``** 

Features currently under development
- Better handling of custom Multiqc v0.7.dev0 via conda to let nextflow deal with conda env creation, and/or
- Migration to Singularity/Docker container
- Complete pipeline with peak calling and IDR analysis.
- Update README with examples & output description
