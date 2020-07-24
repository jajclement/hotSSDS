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
if (params.bamdir){inputType = 'bam'} 
if (params.fqdir){inputType = 'fastq'}

switch (inputType) {
case 'sra':
    Channel
        .fromSRA(params.sra_ids, apiKey:params.ncbi_api_key)
//        .set { sra_ch }
          .set { fq_ch }
    // The process will get the SRA fastq files by pairs or not and will get a small subset
//    process getSRAfiles {
    //    tag { input SRA }
        // this indicates that all output files matching the pattern will be copied into the params.outdir directory.
//        publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.fastq.gz"
        // get input from sra channel. sampleId contains the root names of sample(s) ans reads contains the fastq file(s)
//        input:
//            set val(sampleId), file(reads) from sra_ch
//        output:
//            file '*_1.fastq.gz' into fq1
//            file '*_2.fastq.gz' optional true into fq2        
//        script:
//        """
//            ln -s ${sampleId}_1.fastq.gz ${sampleId}_1.fastq.gz
//            ln -s ${sampleId}_2.fastq.gz ${sampleId}_2.fastq.gz
//        """
//    }   
break
case 'bam':
    Channel
        .fromPath(params.bamdir, checkIfExists:true)
        .set { bam_ch }
    process bam2fastq {
        publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.fastq.gz"
        input:
            file(bam) from bam_ch
        output:
            set val(${bam.baseName}), file('*_1.fastq.gz'), file('*_2.fastq.gz') into fq_ch                                    
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
            samtools sort -n fixmate.bam -o sorted.bam
            bedtools bamtofastq -i sorted.bam -fq ${bam.baseName}_1.fastq -fq2 ${bam.baseName}_2.fastq
            gzip ${bam.baseName}_1.fastq ; gzip ${bam.baseName}_2.fastq
        fi
        """
    }
break
case 'fastq':
    Channel
        .fromFilePairs(params.fqdir, checkIfExists:true)
        .set { fq_ch }
break
}

process fastqc {
publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.fq"
    input:
        set val(sampleId), file(reads) from fq_ch
    output:
        file '*1.fq' into fq1
        file '*2.fq' into fq2
    script:
    """
    zcat ${sampleId}_1.fastq.gz | head -10 > test1.fq
    zcat ${sampleId}_2.fastq.gz | head -10 > test2.fq 
    """
}

