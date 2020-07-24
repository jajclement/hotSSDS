#!/usr/bin/env nextflow
/*
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
//Global parameters
scratch=params.scratch

// Channel to get SRA fastq files from SRA ids given in parameters or config file. Currently no check for availability.
def inputType
if (params.sra_ids){inputType = 'sra'}
//if (params.obj){inputType = 'obj'} 
if (params.bamdir){inputType = 'bam'} 
//if (params.fq1){inputType = 'fastq'}

switch (inputType) {
case 'sra':
    Channel
        .fromSRA(params.sra_ids, apiKey:params.ncbi_api_key)
        .set { sra_ch }

    // The process will get the SRA fastq files by pairs or not and will get a small subset
    process getSRAfiles {
    //    tag { input SRA }
        // this indicates that all output files matching the pattern will be copied into the params.outdir directory.
        publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.fastq.gz"
        // get input from sra channel. sampleId contains the root names of sample(s) ans reads contains the fastq file(s)
        input:
            set val(sampleId), file(reads) from sra_ch
        output:
            file '*1_head.fastq.gz' into fq1
            file '*2_head.fastq.gz' optional true into fq2        
        script:
        """
            zcat ${sampleId}_1.fastq.gz | head -10 > ${sampleId}_1_head.fastq.gz
            zcat ${sampleId}_2.fastq.gz | head -10 > ${sampleId}_2_head.fastq.gz
        """
    }   
break
case 'bam':
    Channel
        .fromPath(params.bamdir, checkIfExists:true)
        .set { bam_ch }
    process bam2fastq {
        publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.fastq.gz"
        input:
            set file(bam) from bam_ch
        output:
            file '*_1.fastq.gz' into fq1                                    
            file '*_2.fastq.gz' optional true into fq2   
        script:
        """
        n=`samtools view -h ${bam} | head -100000 | samtools view -f 1 -S | wc -l`
        if [ \$n -eq 0 ]; then    
            samtools sort -n ${bam} -o sorted.bam
            bedtools bamtofastq -i sorted.bam -fq ${bam.baseName}.fastq
            gzip ${bam.baseName}.fastq
        else
            picard FixMateInformation I=${bam} O=fixmate.bam SORT_ORDER=queryname \
                TMP_DIR=${scratch} VALIDATION_STRINGENCY=LENIENT
            samtools sort -n fixmate -o sorted.bam
            bedtools bamtofastq -i sorted.bam -fq ${bam.baseName}_1.fastq -f2 ${bam.baseName}_2.fastq
            gzip ${bam.baseName}_1.fastq ; gzip ${bam.baseName}_2.fastq
        fi
        """
    }
    break
}











