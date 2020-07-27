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
// Define global variables
def inputType
if (params.sra_ids){inputType = 'sra'}
if (params.bamdir){inputType = 'bam'} 
if (params.fqdir){inputType = 'fastq'}

def outNameStem = "${params.name}.SSDS.${params.genome}"
def tmpNameStem = "${params.name}.tmpFile.${params.genome}"
def generateFastQCScreenConfig_script = "${params.src}/generateFastQCScreenConfig.pl"

// Get input files according to the input types : sra ; bam or fastq
switch (inputType) {
case 'sra':
    Channel
        .fromSRA(params.sra_ids, apiKey:params.ncbi_api_key)
        .set { fq_ch }
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
            set val("${bam.baseName}"), file('*.fastq.gz') into fq_ch                                    
        script:
        """
        picard FixMateInformation I=${bam} O=fixmate.bam SORT_ORDER=queryname \
                    TMP_DIR=${params.scratch} VALIDATION_STRINGENCY=LENIENT
        samtools sort -n fixmate.bam -o sorted.bam
        bedtools bamtofastq -i sorted.bam -fq ${bam.baseName}_1.fastq -fq2 ${bam.baseName}_2.fastq
        gzip ${bam.baseName}_1.fastq ; gzip ${bam.baseName}_2.fastq
        """
    }
break
case 'fastq':
    Channel
        .fromFilePairs(params.fqdir, checkIfExists:true)
        .set { fq_ch }
break
}

process trimFASTQ {
    tag { outNameStem } 
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*_report.txt"
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.html"
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.zip"
    input:
        set val(sampleId), file(reads) from fq_ch
    output:
        set val("${sampleId}_raw"), file(reads) into fqc_ch 
        set val("${sampleId}_trim"),file('*R1.fastq.gz'),file('*R2.fastq.gz') into trim_ch
        file '*_report.txt' into trimmomaticReport
        file '*.html' into trimmedFastqcReport
        file '*.zip' into trimmedFastqcData
    script:
    """
    trimmomatic PE -threads ${task.cpus} ${reads} \
                ${sampleId}_trim_R1.fastq.gz R1_unpaired.fastq.gz \
                ${sampleId}_trim_R2.fastq.gz R2_unpaired.fastq.gz \
                ILLUMINACLIP:${params.adapters}:${params.trim_illuminaclip} SLIDINGWINDOW:${params.trim_slidingwin} MINLEN:${params.trim_minlen} CROP:${params.trim_crop} \
                >& ${sampleId}_trim_${outNameStem}_trimmomatic_report.txt 2>&1
    fastqc -t ${task.cpus} ${sampleId}_trim_R1.fastq.gz ${sampleId}_trim_R2.fastq.gz 
    """
}

process runFASTQC {
    tag { outNameStem }
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.html"
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.zip"
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.png"
    publishDir params.outdir, mode: 'copy', overwrite: false, pattern: "*.txt"
    input:
        set val(sampleId), file(reads) from fqc_ch
    output:
        file '*zip'  into fqcZip
        file '*html' into repHTML
        file '*png'  into repPNG
        file '*txt'  into repTXT
    script:
    """
    fastqc -t ${task.cpus} ${reads}
    perl ${generateFastQCScreenConfig_script} ${params.genomes2screen} ${task.cpus} ${params.genomedir} > fastq_screen.conf 
    fastq_screen --threads ${task.cpus} --force --aligner bwa --conf fastq_screen.conf ${reads}
    """
}

trim_ch
    .splitFastq( by: params.bwaSplitSz, pe:true, file:true)
    .set { splitFQ_ch }

process bwaAlign {
    tag { outNameStem }
    input:
        set val(sampleId), file(fqR1), file(fqR2) from splitFQ_ch
    output:
        file 'bar*bam' into multiBAMaln, multiBAM2merge
    script:
    def nm = new Random().with {(1..30).collect {(('a'..'z')).join()[ nextInt((('a'..'z')).join().length())]}.join()} 
    def bamOut = 'bar_' + nm + '.bam'
    """
    if [ ${params.trim_cropR1} != ${params.trim_cropR2} ]
    then
        fastx_trimmer -f 1 -l ${params.trim_cropR1} -i ${fqR1} -o ${tmpNameStem}.R1.fastq
        fastx_trimmer -f 1 -l ${params.trim_cropR2} -i ${fqR2} -o ${tmpNameStem}.R2.fastq
    else
        mv ${fqR1} ${tmpNameStem}.R1.fastq
        mv ${fqR2} ${tmpNameStem}.R2.fastq
    fi
    ${params.custom_bwa} aln -t ${task.cpus} ${params.genome_fasta} ${tmpNameStem}.R1.fastq \
            > ${tmpNameStem}.R1.sai
    ${params.custom_bwa_ra} aln -t ${task.cpus} ${params.genome_fasta} ${tmpNameStem}.R2.fastq \
            > ${tmpNameStem}.R2.sai
    ${params.custom_bwa} sampe ${params.genome_fasta} ${tmpNameStem}.R1.sai ${tmpNameStem}.R2.sai ${tmpNameStem}.R1.fastq \
            ${tmpNameStem}.R2.fastq >${tmpNameStem}.unsorted.sam
    picard SamFormatConverter I=${tmpNameStem}.unsorted.sam O=${tmpNameStem}.unsorted.tmpbam VALIDATION_STRINGENCY=LENIENT
    picard SortSam I=${tmpNameStem}.unsorted.tmpbam O=${bamOut} SO=coordinate VALIDATION_STRINGENCY=LENIENT
    samtools index ${bamOut}
    ls -l >list.tab
    """
}


