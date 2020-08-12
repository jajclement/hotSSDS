## **SSDS nextflow pipeline version 1.8_NF_pa**
### **Welcome**
**This is [SSDS pipeline by Kevin Brick](https://github.com/kevbrick/SSDSnextflowPipeline) updated and adapted to IGH cluster.**
See [initial paper](https://genome.cshlp.org/content/22/5/957.long) and [technical paper](https://www.sciencedirect.com/science/article/pii/S0076687917303750?via%3Dihub).
The pipeline uses [Nextflow]( https://www.nextflow.io/).
Briefly, the update from SSDS pipeline version 1.8_NF included **conda profile**, input and trimming process modification, and global nextflow homogeneisation.

*Please report any bug occured during the installation to pauline.auffret at igh dot cnrs dot fr*


## **How to run the pipeline on IGH cluster**
### Requirements
Get [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) if not installed on your system.
Then install **Nextflow version 20.04.1** through conda :
````sh
conda create -name nextflow-dev
conda activate nextflow-dev
conda install -c bioconda nextflow=20.04.1
conda deactivate 
````
You may need to create a sbatch script to run these commands
### 1. Get the pipeline
First you need to download the pipeline in your working directory
````sh
cd /home/${USER}/work
git clone https://gitlab.igh.cnrs.fr/pauline.auffret/SSDSnextflowPipeline.git
cd SSDSnextflowPipeline
````
### 2. Pipeline configuration 
There are currently 2 configuration files :
- ````./conf/igh.config```` contains cluster requirements. You don't need to edit this file unless you want to custom the requirements for CPU/memory usage and compute queue (see [IGH cluster documentation](https://kojiki.igh.cnrs.fr/doku.php?id=cluster,))
- ````./nextflow.config```` contains default pipeline parameters. You *can* edit this file but it is recommendend to use ````run_pipeline.sh```` script to use your own parameters (this will overwrite the default configuration), see section 4.

### 3. Input data
The pipeline will only process **paired-end data**.
Currently 3 input data formats are supported : 
- fastq(.gz) ````--fqdir input_data/raw_data/*{R1,R2}.fastq.gz````
- bam ````--bamdir input_data/raw_data/*.bam````
- SRA identifiers ````--sra_ids=['ERR908507', 'ERR908506']````
The reference genomes should be contained in a folder following a very specific structure (you can use soft links)

### 4. Run the pipeline !


### Test data
You may want to test the installtion before going with your own date. If so, just...


Features currently under development
- Access to genomes nf-core style
- Migration to Singularity/socker containeur
