#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                        SSDS Pipeline version 1.8_NF
			Author : Kevin Brick
========================================================================================
 SSDS nextflow pipeline
 #### Homepage / Documentation
 https://github.com/kevbrick/SSDSnextflowPipeline
 Adapted from version 1.8_NF (Pauline Auffret, 2020)
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

	nextflow run ...
	or use the bash pipeline

    Mandatory arguments:
		TO DO.
"""
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                SET UP CONFIGURATION VARIABLES                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}
SCRATCH=params.scratch
fastqdir=params.fastqdir
bamdir=params.bamdir

// Get the input type : sra ID ; bam files ; fastq PE files ; fastq SR files
def inputType 
    if (params.sra){inputType = 'sra'}                   
    if (params.bamdir){inputType = 'bam'}
    if (params.fastqdir){inputType = 'fq'}

switch (inputType) {
    case 'bam':
        Channel
            .fromPath(params.bamdir, checkIfExists:true) 
            set { bam_ch }

        process bamToFastq {
            tag { bam } 
            publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*expttype"
            input:
               file bam from bam_ch

            output:
               fastq_ch

            script:
               """
               n=`samtools view -h ${bam} | head -100000 | samtools view -f 1 -S | wc -l`
               if [ \$n -eq 0 ]; then
                   echo 1 >isSR.expttype
                   picard SortSam I=${bam} O=querynameSort.bam SO=queryname TMP_DIR=$SCRATCH VALIDATION_STRINGENCY=LENIENT 2>b2fq.err
                   



    case 'sra':
        Channel                                                                                                   .fromSRA(params.sra, checkIfExists:true)                                                              .set { sra_ch }
        process getSRAfiles {
            tag { params.sra }
            input:
                set sampleId, file(reads) from sra_ch
            output:
                  publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*expttype"

}
}
process input_sra {
    output:
    
Channel
    .fromSRA(params.sra, checkIfExists:true)
    .set { sra_ch }

Channel
    .fromPath(params.bamdir, checkIfExists:true)
    set { bam_ch }

Channel
    .fromFilePairs("$fastqdir/*_{R1,R2}*.fastq.gz", checkIfExists:true)
    .set { fastqPE_ch }

process bam2fq {
    input:
    bam_ch = Channel.fromPath("$bamdir/*.bam", checkIfExists:true)

process foo {
  input:
  set sampleId, file(reads) from samples_ch

  script:
  """
  echo your_command --sample $sampleId --reads $reads
  """
}

