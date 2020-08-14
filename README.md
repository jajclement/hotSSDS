## **SSDS nextflow pipeline version 1.8_NF_pa**
**/!\ Work in progress /!\**
### **Welcome**
**This is [SSDS pipeline by Kevin Brick](https://github.com/kevbrick/SSDSnextflowPipeline) updated and adapted to IGH cluster.**
See [initial paper](https://genome.cshlp.org/content/22/5/957.long) and [technical paper](https://www.sciencedirect.com/science/article/pii/S0076687917303750?via%3Dihub).
The pipeline uses [Nextflow]( https://www.nextflow.io/).
Briefly, the update from SSDS pipeline version 1.8_NF included **conda profile**, input and trimming process modification, and global nextflow homogeneisation.

*Please report any bug occured during the installation to pauline.auffret at igh dot cnrs dot fr*


## **How to run the pipeline on IGH cluster**
### Requirements
Get [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) if not installed on your system.

### 1. Get the pipeline and set up conda environment
First you need to download the pipeline in your working directory
````sh
cd /home/${USER}/work
git clone https://gitlab.igh.cnrs.fr/pauline.auffret/SSDSnextflowPipeline.git
cd SSDSnextflowPipeline
````
Then install pipeline conda environment through sbatch script :
```` 
sbatch -p computepart -J "install_conda_env" --export=ALL --mem 5G -t 5-0:0 --mem-per-cpu=1000 --wrap "bash src/install_pipeline.sh"
```` 
This will create 2 conda environments : nextflow-dev and SSDSnextflowPipeline.
If you have existing conda environments with the same names, you need to edit the install_pipeline.sh
*This is a temporary setting*

### 2. Pipeline configuration 
There are currently 2 configuration files :
- ````./conf/igh.config```` contains cluster resources requirements & reference genomes info. You don't need to edit this file unless you want to custom the requirements for CPU/memory usage and compute queue (see [IGH cluster documentation](https://kojiki.igh.cnrs.fr/doku.php?id=cluster,))
- ````./nextflow.config```` contains default pipeline parameters. You *can* edit this file but it is recommendend to use ````run_pipeline.sh```` script to use your own parameters (this will overwrite the default configuration), see section 4.

### 3. Input data
The pipeline will only process **paired-end data**.
Currently 3 input data formats are supported : 
- fastq(.gz) ````--fqdir input_data/raw_data/*{R1,R2}.fastq.gz````
- bam ````--bamdir input_data/raw_data/*.bam````
- SRA identifiers ````--sra_ids=['ERR908507', 'ERR908506']````

The reference genome should be in the ``/poolzvs/genomes`` directory on IGH cluster. Currently, the available genomes are mm10, hg19, hg38, sacCer2, sacCer3, dm3, dm6.
If so, you can use ````-- genome 'mm10'```` et that's it for the genome parameters.
If you want to use another reference, you will need to set : **not tested**
- path to genome ````--genomedir /path/to/genome````
- path to fasta file. Indexes for BWA SHOULD EXIST in the same directory ````--genome_fasta /path/to/genome.fa````
- the name of the genome ````--genome_name mm11````
- the path to the fai index ````--fai /path/to/genome.fai.fai````

### 4. Run the pipeline !
To see all the parameters, please run :
````
cd /home/${USER}/work/SSDSnextflowPipeline
conda activate nextflow-dev
nextflow run main.nf --help
````
You can either run the pipeline directly through the command line :
````
cd /home/${USER}/work/SSDSnextflowPipeline
conda activate nextflow-dev
sbatch -p computepart -J SSDSnextflowPiepline --export=ALL --mem 5G -t 5-0:0 --mem-per-cpu=1000 \
	--wrap "nextflow run main.nf -c conf/igh.config --fqdir path/to/your/data/*{R1,R2}.fastq.gz --name your_analysis_name --genome mm10"
````
or use ``run_pipeline.sh`` script : in this case make sure to edit with your parameters, then
````
bash run_pipeline.sh
````

### Test data
You may want to test the installation before going with your own date. If so, just use ````--fqdir tests/fastq/*{R1,R2}.fastq"````
e.g.
````
cd /home/${USER}/work/SSDSnextflowPipeline
conda activate nextflow-dev
sbatch -p computepart -J SSDSnextflowPiepline --export=ALL --mem 5G -t 5-0:0 --mem-per-cpu=1000 \
	--wrap "nextflow run main.nf -c conf/igh.config --fqdir tests/fastq/*{R1,R2}.fastq --name your_analysis_name --genome mm10"
````
or use ``run_pipeline.sh`` script (make sur to edit it with the relevant paramters before use) :
````
bash run_pipeline.sh
````

Features currently under development
- Better handling of custom Multiqc v0.7.dev0 via conda to let nextflow deal with conda env creation, and/or
- Migration to Singularity/Docker container
- Complete pipeline with peak calling and IDR analysis.
- Update README with examples & output description
