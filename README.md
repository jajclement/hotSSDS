# SSDS nextflow pipeline : 
Initial paper: Khil et al. Genome Research 2012

Technical paper: Brick, Pratto et al., Methods in Enzymology 2018

This nextflow pipeline is configured to work on a SLURM-based HPC with modules. It can be relatively easily configured to run on other systems (see nextflow documentation : https://www.nextflow.io/). 

# Requirements:
Tools / programs :
- bedtools	2.25.0
- bwa	0.7.17
- deeptools	3.0.1
- fastqc	0.11.8
- fastqtools	0.8
- fastxtoolkit	0.0.14
- java	1.8.0_92
- nextflow	0.30.2
- picard	2.9.2
- python	2.7
- samtools	1.8
- sratoolkit	2.9.2
- trimgalore	0.4.5
- ucsc	365

# Global variables required: 
$NXF_PIPEDIR   : Path to folder containing SSDSPipeline_1.6.groovy

$NXF_GENOMES   : Path to folder containing reference genomes for alignment
                 ** This folder requires a very specific structure (see below) **

$SLURM_JOBID   : Specifies the temporary subfolder to use 

### NXF_GENOMES Folder structure
Each reference genome should be contained in a separate folder (i.e. $NXF_GENOMES/mouse_mm10). The sub-structure within this folder should be as follows: 

$NXF_GENOMES/<genome>/genome.fa                : Genome fasta file

$NXF_GENOMES/<genome>/genome.fa.fai            : Index of genome fasta file (samtools faidx)

$NXF_GENOMES/<genome>/genome.dict              : Sequence dictionary for genome fasta file (use picard CreateSequenceDictionary)

$NXF_GENOMES/<genome>/BWAIndex/version0.7.10/  : BWA 0.7 index files (should also contain soft links to the three files above)



# Other requirements
The pipeline requires a high-level temporary folder called /lscratch. On a SLURM-based HPC, each job is assigned a global id ($SLURM_JOBID) and this is appended to the temp folder name for each process. This is currently hard-coded. Thus, there is a requirement for :

/lscratch folder for temporary files
SLURM_JOBID global variable for each HPC job.

Basic execution is as follows: 

## RUN ON LOCAL MACHINE
nextflow run -c $NXF_PIPEDIR/nextflow.local.config $NXF_PIPEDIR/SSDSPipeline_1.6.groovy \
    --fq1 $NXF_PIPEDIR/tests/fastq/ssdsLong.100k.R1.fastq \
    --fq2 $NXF_PIPEDIR/tests/fastq/ssdsLong.100k.R2.fastq \
    --r1Len 36 \
    --r2Len 40 \
    --genome mm10 \
    --outName testSSDS1.6 \
    --outdir SSDS1.6_test \
    -with-trace -with-timeline

## RUN ON SLURM CLUSTER
nextflow run -c $NXF_PIPEDIR/nextflow.config $NXF_PIPEDIR/SSDSPipeline_1.6.groovy \
    --fq1 $NXF_PIPEDIR/tests/fastq/ssdsLong.100k.R1.fastq \
    --fq2 $NXF_PIPEDIR/tests/fastq/ssdsLong.100k.R2.fastq \
    --r1Len 36 \
    --r2Len 40 \
    --genome mm10 \
    --outName testSSDS1.6 \
    --outdir SSDS1.6_test \
    -with-trace -with-timeline

The tests folder contains small fastq files from an SSDS experiment in mouse. The test should take a few minutes on a local machine (64Gb RAM, 16 Cores).




